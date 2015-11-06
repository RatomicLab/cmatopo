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

using namespace cma;
using namespace std;

const double DEFAULT_TOLERANCE = 1.;

namespace cma {
    GEOSContextHandle_t hdl;
}

int main(int argc, char **argv)
{
    initGEOS(geos_message_function, geos_message_function);
    OGRRegisterAll();

    GEOSHelper geos = GEOSHelper();
    assert (hdl != NULL);

    GEOSWKBReader* wkbr = geos.reader();
    assert (wkbr != NULL);

    // vm with Quebec only: postgresql://postgres@192.168.56.101/postgis
    // vm with the world: postgresql://pgsql@pg/cmdb
    PG db("postgresql://postgres@localhost/postgis");
    if (!db.connected()) {
        cerr << "Could not connect to PostgreSQL." << endl;
        return 1;
    }

    GEOSGeometry* world_extent = world_geom();

    std::vector<zoneInfo*> zones;
    prepare_zones(db, geos, world_extent, zones, 10);
    GEOSGeom_destroy_r(hdl, world_extent);
    // write_zones("output/test.shp", zones, true);

    sort(zones.begin(), zones.end(), [](zoneInfo* a, zoneInfo* b) {
        return a->second > b->second;
    });

    int zc = 0;
    for (zoneInfo* zone : zones) {
        if (zone->second == 0) {
            continue;
        }

        ++zc;

        if (zone->second > 5000) continue;

        chrono::time_point<chrono::system_clock> start, end;
        start = chrono::system_clock::now();

        Topology* topology = new Topology(geos);

        linesV lines;
        if (!db.get_lines_within(zone->first, lines)) {
            assert (false);
        }

        assert (zone->second == lines.size());
        cout << lines.size() << endl;

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

        cout << "finished computation at " << std::ctime(&end_time)
             << "elapsed time: " << elapsed_seconds.count() << "s\n";

        delete topology;
        break;
    }

    finishGEOS();

    return 0;

    int count = 0;
    PGresult* res = db.query("SELECT line from way;", true);
    while (res != NULL && !PQgetisnull(res, 0, 0)) {
        ++count;

        char* line = PQgetvalue(res, 0, 0);

        GEOSGeometry* geom = GEOSWKBReader_readHEX_r(hdl, wkbr, (unsigned char*)line, strlen(line));
        if (!geom) {
            cerr << "Invalid geometry received." << endl;
            return 1;
        }

        if (GEOSGeomTypeId_r(hdl, geom) != GEOS_LINESTRING) {
            char* geom_type = GEOSGeomType_r(hdl, geom);
            cerr << "Skipping invalid geometry of type " << geom_type << endl;
            GEOSFree_r(hdl, geom_type);
        }

        GEOSGeom_destroy_r(hdl, geom);

        PQclear(res);
        res = db.next_result();
    }

    cout << count << endl;
}
