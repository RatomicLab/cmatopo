// icpc -g -openmp -O3 -std=c++11 -I.. geosidx.cpp ../utils.o -o geosidx -lgeos_c && DYLD_LIBRARY_PATH=/opt/intel/compilers_and_libraries/mac/lib ./geosidx

#include <cassert>
#include <iostream>
#include <geos_c.h>

#include <utils.h>

using namespace cma;
using namespace std;

namespace cma {
    GEOSContextHandle_t hdl;
}

void callback(void *item, void *userdata)
{
    GEOSGeometry* geom_contaning = (GEOSGeometry*)item;
    GEOSGeometry* geom_contained = (GEOSGeometry*)userdata;

    GEOSWKTWriter* w = GEOSWKTWriter_create_r(hdl);
    char* s1 = GEOSWKTWriter_write_r(hdl, w, geom_contaning);
    char* s2 = GEOSWKTWriter_write_r(hdl, w, geom_contained);

    cout << s1 << " intersects with " << s2 << endl;

    GEOSFree_r(hdl, s1);
    GEOSFree_r(hdl, s2);
    GEOSWKTWriter_destroy_r(hdl, w);
}

int main(int /*argc*/, char** /*argv*/)
{
    GEOSHelper geos = GEOSHelper();
    assert (hdl != NULL);

    GEOSWKTReader* wktr = geos.text_reader();
    assert (wktr != NULL);

    GEOSSTRtree* index = GEOSSTRtree_create_r(hdl, 10);
    assert (index != NULL);

    GEOSGeometry* bigbox = GEOSWKTReader_read_r(hdl, wktr, "LINESTRING(0 0, 10 0, 10 10, 0 10, 0 0)");
    GEOSSTRtree_insert_r(hdl, index, bigbox, (void*)bigbox);

    GEOSGeometry* contained_box = GEOSWKTReader_read_r(hdl, wktr, "LINESTRING(2 2, 8 2, 8 8, 2 8, 2 2)");
    GEOSSTRtree_insert_r(hdl, index, contained_box, contained_box);

    GEOSGeometry* intersect_box = GEOSWKTReader_read_r(hdl, wktr, "LINESTRING(-5 5, 5 5, 5 15, -5 15, -5 5)");
    GEOSSTRtree_insert_r(hdl, index, intersect_box, intersect_box);

    GEOSGeometry* outside_box = GEOSWKTReader_read_r(hdl, wktr, "LINESTRING(-100 -100, -50 -100, -50 -50, -100 -50, -100 -100)");
    GEOSSTRtree_insert_r(hdl, index, outside_box, outside_box);

    GEOSSTRtree_query_r(hdl, index, contained_box, &callback, contained_box);
    GEOSSTRtree_query_r(hdl, index, outside_box, &callback, outside_box);

    GEOSSTRtree_destroy_r(hdl, index);
    GEOSGeom_destroy_r(hdl, bigbox);
    GEOSGeom_destroy_r(hdl, contained_box);
    return 0;
}
