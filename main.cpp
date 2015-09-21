#include <cassert>
#include <iostream>

#include <geos_c.h>
#include <ogrsf_frmts.h>

#include <pg.h>
#include <utils.h>
#include <zones.h>

using namespace cma;
using namespace std;

namespace cma {
    GEOSContextHandle_t hdl;
}

int main(int argc, char **argv)
{
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
    write_zones("output/test.shp", zones, true);

    return 0;
}
