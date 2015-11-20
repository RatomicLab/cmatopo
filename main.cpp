#include <ctime>
#include <chrono>
#include <cassert>
#include <iostream>

#include <geos_c.h>
#include <ogrsf_frmts.h>

#include <pg.h>
#include <utils.h>
#include <zones.h>
#include <topology.h>

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>

using namespace cma;
using namespace std;

using namespace boost::mpi;

namespace cma {
    GEOSContextHandle_t hdl;
}

int main(int argc, char **argv)
{
    environment env;
    communicator world;

    initGEOS(geos_message_function, geos_message_function);
    OGRRegisterAll();

    unique_ptr<GEOSHelper> geos(new GEOSHelper());
    assert (hdl != NULL);

    // vm with Quebec only: postgresql://postgres@192.168.56.101/postgis
    // vm with the world: postgresql://pgsql@pg/cmdb
    // PG db("postgresql://postgres@localhost/postgis");
    PG db("postgresql://laurent@localhost/cmatopo");
    if (!db.connected()) {
        cerr << "Could not connect to PostgreSQL." << endl;
        return 1;
    }

    int line_count = db.get_line_count();
    assert (line_count >= 0);

    vector<int>* zonesPerProcess = new vector<int>(world.size());
    vector< list< pair<int, string> > >* processZones =
        new vector< list< pair<int, string> > >(world.size());

    assert (zonesPerProcess->size() == processZones->size());

    for (int i = 0; i < zonesPerProcess->size(); ++i) {
        (*zonesPerProcess)[i] = 0;
    }

    int processingLineCount = 0;
    if (world.rank() == 0) {
        std::vector<zoneInfo*> zones;

        GEOSGeometry* world_extent = world_geom();
        cout << "count: " << prepare_zones(*geos, world_extent, zones, 10) << endl;
        GEOSGeom_destroy_r(hdl, world_extent);
        // write_zones("output/test.shp", zones, true);

        sort(zones.begin(), zones.end(), [](zoneInfo* a, zoneInfo* b) {
            return a->second > b->second;
        });

        int zoneId = 0;
        for (const zoneInfo* z : zones) {
            int minProcess = distance(
                zonesPerProcess->begin(),
                min_element(
                    zonesPerProcess->begin(),
                    zonesPerProcess->end()
                )
            );
            assert (minProcess < zonesPerProcess->size());
            (*zonesPerProcess)[minProcess] += z->second;

            GEOSGeometry* zoneGeom = OGREnvelope2GEOSGeom(z->first);
            (*processZones)[minProcess].push_back(
                make_pair(zoneId++, geos->as_string(zoneGeom))
            );
            GEOSGeom_destroy_r(hdl, zoneGeom);
            delete z;
        }
        zones.clear();

        for (int i = 0; i < zonesPerProcess->size(); ++i) {
            cout << "[0] Process " << i << " will process " << (*zonesPerProcess)[i]
                 << " lines." << endl;
            processingLineCount += (*zonesPerProcess)[i];
        }

        cout << "Will process " << processingLineCount << ", leaving " << line_count-processingLineCount
             << " orphans (" << (line_count-processingLineCount)/float(processingLineCount)*100 << "%)" << endl;
    }

    list< pair<int, string> > myZones;
    scatter(world, *processZones, myZones, 0);

    processZones->clear();
    zonesPerProcess->clear();

    delete processZones;
    delete zonesPerProcess;

    vector<Topology*> myTopologies;
    unique_ptr< set<int> > myOrphans(new set<int>());

    for (auto& zone : myZones) {
        int zoneId = zone.first;
        const string& hexWKT = zone.second;

        chrono::time_point<chrono::system_clock> start, end;
        start = chrono::system_clock::now();

        GEOSGeometry* zoneGeom =
            GEOSWKTReader_read_r(hdl, geos->text_reader(), hexWKT.c_str());

        linesV lines;
        if (!db.get_lines(zoneGeom, lines, true)) {
            assert (false);
        }

        if (!db.get_line_ids(zoneGeom, *myOrphans, false)) {
            assert (false);
        }

        GEOSGeom_destroy_r(hdl, zoneGeom);

        cout << "[" << world.rank() << "] processing zone #" <<  zoneId
             << " (" << lines.size() << " lines, "
             << myOrphans->size() << " orphans)" << endl;

        if (lines.size() == 0) {
            continue;
        }

        Topology* topology = new Topology(*geos);

        int lc = 0;
        for (GEOSGeometry* line : lines) {
            vector<int> edgeIds;
            try {
                topology->TopoGeo_AddLineString(line, edgeIds, DEFAULT_TOLERANCE);
                topology->commit();
            }
            catch (const invalid_argument& ex) {
                cerr << geos->as_string(line) << ": " << ex.what() << endl;
                topology->rollback();
            }

            GEOSGeom_destroy_r(hdl, line);

            if (++lc % 100 == 0) {
                cout << lc << endl;
            }
        }
        lines.clear();

        end = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds = end-start;
        time_t end_time = chrono::system_clock::to_time_t(end);

        cout << "[" << world.rank() << "] finished computation of zone #" << zoneId
             << " at " << std::ctime(&end_time) << ","
             << " elapsed time: " << elapsed_seconds.count() << "s" << endl;

        myTopologies.push_back(topology);
        myZones.pop_front();
    }

    Topology* mainTopology = myTopologies[0];

    if (world.rank() == 0) {
        // merge my own topologies first
        for (Topology* otherTopology : myTopologies) {
            if (mainTopology == otherTopology) {
                continue;
            }
            merge_topologies(*mainTopology, *otherTopology);
        }

        // TODO: gather all topologies (one process at a time?)
        /*
        for (int rank = 1; rank < world.rank(); ++rank) {
            vector<Topology> topologies;
            world.recv(rank, 0, topologies);
        }
        */

        mainTopology->rebuild_indexes();

        unique_ptr< vector< set<int> > > allOrphansV(new vector< set<int> >());
        gather(world, *myOrphans, *allOrphansV, 0);

        unique_ptr< set<int> > allOrphans(new set<int>());
        for (auto& otherOrphans : *allOrphansV) {
            allOrphans->insert(otherOrphans.begin(), otherOrphans.end());
        }

        cout << "Found a total of " << allOrphans->size() << " orphans." << endl;

        vector<int> edgeIds;
        for (int _id : *allOrphans) {
            GEOSGeometry* line = db.get_line(_id);
            cout << _id << ",";
            mainTopology->TopoGeo_AddLineString(line, edgeIds, DEFAULT_TOLERANCE);
        }

        // FIXME: assert (allOrphans->size() == line_count-processingLineCount);

        mainTopology->output();
    }
    else {
        // TODO: world.send to rank 0

        gather(world, *myOrphans, 0);
    }

    for (Topology* t : myTopologies) {
        delete t;
    }
    myTopologies.clear();

    finishGEOS();

    return 0;
}
