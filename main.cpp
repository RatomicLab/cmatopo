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
    // GEOSGeometry* env = GEOSEnvelope_r(hdl, world_extent);
    // cout << env << endl;
    // env = GEOSEnvelope_r(hdl, world_extent);
    // cout << env << endl << GEOSEnvelope_r(hdl, world_extent) << endl;
    //
    // return 0;

    std::vector<zoneInfo*> zones;
    prepare_zones(db, geos, world_extent, zones, 10);
    GEOSGeom_destroy_r(hdl, world_extent);
    // write_zones("output/test.shp", zones, true);

    Topology* topology = new Topology(geos);

    for (zoneInfo* zone : zones) {
        if (zone->second == 0) {
            continue;
        }

        linesV lines;
        if (!db.get_lines_within(zone->first, lines)) {
            assert (false);
        }

        assert (lines.size() > 0);

        for (GEOSGeometry* line : lines) {
            vector<int> edgeIds;
            try {
                topology->TopoGeo_AddLineString(line, edgeIds, DEFAULT_TOLERANCE);
            }
            catch (const invalid_argument& ex) {
                cerr << geos.as_string(line) << ": " << ex.what() << endl;
            }
        }

        topology->output_nodes();
        topology->output_edges();
    }

    delete topology;

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
