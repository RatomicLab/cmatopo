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
#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>

using namespace cma;
using namespace std;

using namespace boost::mpi;

const double DEFAULT_TOLERANCE = 1.;

namespace cma {
    GEOSContextHandle_t hdl;
}

int main(int argc, char **argv)
{
    environment env;
    communicator world;

    initGEOS(geos_message_function, geos_message_function);
    OGRRegisterAll();

    GEOSHelper geos = GEOSHelper();
    assert (hdl != NULL);

    // vm with Quebec only: postgresql://postgres@192.168.56.101/postgis
    // vm with the world: postgresql://pgsql@pg/cmdb
    PG db("postgresql://postgres@localhost/postgis");
    if (!db.connected()) {
        cerr << "Could not connect to PostgreSQL." << endl;
        return 1;
    }

    vector<int>* zonesPerProcess = new vector<int>(world.size());
    vector< list< pair<int, string> > >* processZones =
        new vector< list< pair<int, string> > >(world.size());

    assert (zonesPerProcess->size() == processZones->size());

    for (int i = 0; i < zonesPerProcess->size(); ++i) {
        (*zonesPerProcess)[i] = 0;
    }

    if (world.rank() == 0) {
        std::vector<zoneInfo*> zones;

        GEOSGeometry* world_extent = world_geom();
        prepare_zones(geos, world_extent, zones, 10);
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
                make_pair(zoneId++, geos.as_string(zoneGeom))
            );
            GEOSGeom_destroy_r(hdl, zoneGeom);
        }

        for (int i = 0; i < zonesPerProcess->size(); ++i) {
            cout << "[0] Process " << i << " will process " << (*zonesPerProcess)[i]
                 << " lines." << endl;
        }
    }

    list< pair<int, string> > myZones;
    scatter(world, *processZones, myZones, 0);

    processZones->clear();
    zonesPerProcess->clear();

    delete processZones;
    delete zonesPerProcess;

    for (auto& zone : myZones) {
        int zoneId = zone.first;
        const string& hexWKT = zone.second;

        chrono::time_point<chrono::system_clock> start, end;
        start = chrono::system_clock::now();

        GEOSGeometry* zoneGeom =
            GEOSWKTReader_read_r(hdl, geos.text_reader(), hexWKT.c_str());

        linesV lines;
        if (!db.get_lines_within(zoneGeom, lines)) {
            assert (false);
        }

        cout << "[" << world.rank() << "] processing zone #" <<  zoneId
             << " (" << lines.size() << " lines)";

        if (lines.size() == 0) {
            continue;
        }

        Topology* topology = new Topology(geos);

        int lc = 0;
        for (GEOSGeometry* line : lines) {
            vector<int> edgeIds;
            try {
                topology->TopoGeo_AddLineString(line, edgeIds, DEFAULT_TOLERANCE);
            }
            catch (const invalid_argument& ex) {
                cerr << geos.as_string(line) << ": " << ex.what() << endl;
            }

            if (++lc % 100 == 0) {
                cout << lc << endl;
            }
        }

        end = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds = end-start;
        time_t end_time = chrono::system_clock::to_time_t(end);

        topology->output_nodes();
        topology->output_edges();
        topology->output_faces();

        cout << "[" << world.rank() << "] finished computation of zone #" << zoneId
             << " at " << std::ctime(&end_time) << ","
             << " elapsed time: " << elapsed_seconds.count() << "s" << endl;

        delete topology;
    }

    finishGEOS();

    return 0;
}
