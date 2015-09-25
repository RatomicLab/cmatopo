#include <st.h>

#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>

extern "C" {
    #include <liblwgeom.h>
    #include <lwgeom_geos.h>
}

using namespace std;

namespace cma {

const int DIST_MIN = 1;

bool ST_Equals(const GEOSGeom g1, const GEOSGeom g2)
{
    return GEOSWithin(g1, g2) == 1 && GEOSWithin(g2, g1) == 1;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_dwithin(PG_FUNCTION_ARGS)
 */
bool ST_DWithin(const GEOSGeom g1, const GEOSGeom g2, double tolerance)
{
    assert (tolerance >= 0.);

    LWGEOM* lwgeom1 = GEOS2LWGEOM(g1, 0);
    LWGEOM* lwgeom2 = GEOS2LWGEOM(g2, 0);

    double distance = lwgeom_mindistance2d_tolerance(lwgeom1, lwgeom2, tolerance);

    lwgeom_free(lwgeom1);
    lwgeom_free(lwgeom2);

    return distance < tolerance;
}

bool ST_IsEmpty(const GEOSGeom geom)
{
    return GEOSisEmpty_r(hdl, geom) == 1;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_closestpoint(PG_FUNCTION_ARGS)
 */
bool ST_Contains(const GEOSGeom g1, const GEOSGeom g2)
{
    return GEOSContains(g1, g2);
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_same(PG_FUNCTION_ARGS)
 */
bool ST_OrderingEquals(const GEOSGeom g1, const GEOSGeom g2)
{
    if (GEOSGeomTypeId(g1) != GEOSGeomTypeId(g2)) {
        return false;
    }

    LWGEOM* lwg1 = GEOS2LWGEOM(g1, 0);
    LWGEOM* lwg2 = GEOS2LWGEOM(g2, 0);

    bool ret = lwgeom_same(lwg1, lwg2);

    lwgeom_free(lwg1);
    lwgeom_free(lwg2);

    return ret;
}

double ST_X(const GEOSGeom geom)
{
    double x;
    GEOSGeomGetX(geom, &x);
    return x;
}

double ST_Y(const GEOSGeom geom)
{
    double y;
    GEOSGeomGetY(geom, &y);
    return y;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_azimuth(PG_FUNCTION_ARGS)
 */
double ST_Azimuth(const GEOSGeometry* g1, const GEOSGeometry* g2)
{
    assert (g1 && GEOSGeomTypeId_r(hdl, g1) == GEOS_POINT);
    assert (g2 && GEOSGeomTypeId_r(hdl, g2) == GEOS_POINT);

    LWGEOM* lwgeom1 = GEOS2LWGEOM(g1, 0);
    LWGEOM* lwgeom2 = GEOS2LWGEOM(g2, 0);

    LWPOINT* lwpt1 = lwgeom_as_lwpoint(lwgeom1);
    LWPOINT* lwpt2 = lwgeom_as_lwpoint(lwgeom2);

    assert (lwpt1->srid == lwpt2->srid);

    POINT2D p1, p2;
    if (!getPoint2d_p(lwpt1->point, 0, &p1) || !getPoint2d_p(lwpt2->point, 0, &p2)) {
        lwgeom_free(lwgeom1);
        lwgeom_free(lwgeom2);
        return -1.;
    }

    double az;
    if (!azimuth_pt_pt(&p1, &p2, &az)) {
        return -1.;
    }

    lwgeom_free(lwgeom1);
    lwgeom_free(lwgeom2);

    return az;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_mindistance2d(PG_FUNCTION_ARGS)
 */
double ST_Distance(const GEOSGeom g1, const GEOSGeom g2)
{
    LWGEOM* lwgeom1 = GEOS2LWGEOM(g1, 0);
    LWGEOM* lwgeom2 = GEOS2LWGEOM(g2, 0);

    double distance = lwgeom_mindistance2d(lwgeom1, lwgeom2);

    lwgeom_free(lwgeom1);
    lwgeom_free(lwgeom2);

    return distance;
}

double _ST_MinTolerance(const GEOSGeom geom)
{
    vector<double> bbox;
    bounding_box(geom, bbox);

    int _max = (int)*max_element(bbox.begin(), bbox.end());
    int val = _max != 0 ? _max : 1;

    return 3.6 * pow(10,  -(15 - log(val)));
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum ST_Snap(PG_FUNCTION_ARGS)
 */
GEOSGeom ST_Snap(const GEOSGeom g1, const GEOSGeom g2, double tolerance)
{
    assert (g1 && g2);

    int srid = GEOSGetSRID_r(hdl, g1);

    GEOSGeometry* g3 = GEOSSnap_r(hdl, g1, g2, tolerance);

    if (g3) {
        GEOSSetSRID_r(hdl, g3, srid);
    }

    return g3;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_geos.c):
 *   Datum ST_Split(PG_FUNCTION_ARGS)
 */
GEOSGeom ST_Split(const GEOSGeometry* in, const GEOSGeometry* blade_in)
{
    LWGEOM* lwgeom_in  = GEOS2LWGEOM(in, 0);
    LWGEOM* lwblade_in = GEOS2LWGEOM(blade_in, 0);
    LWGEOM* lwgeom_out;

    lwgeom_out = lwgeom_split(lwgeom_in, lwblade_in);
    lwgeom_free(lwgeom_in);
    lwgeom_free(lwblade_in);

    GEOSGeom ret = LWGEOM2GEOS(lwgeom_out);
    lwgeom_free(lwgeom_out);

    return ret;
}

GEOSGeom ST_PointN(const GEOSGeometry* line, int index)
{
    assert (line && GEOSGeomTypeId_r(hdl, line) == GEOS_LINESTRING);
    assert (index >= 0 && index < GEOSGeomGetNumPoints_r(hdl, line));
    return GEOSGeomGetPointN_r(hdl, line, index);
}

GEOSGeom ST_Collect(GEOSGeometry* g1, GEOSGeometry* g2)
{
    GEOSGeometry* geoms[] = {g1, g2};

    // After this point, g1 and g2 should not be freed.
    // Only the created collection should.

    return GEOSGeom_createCollection_r(hdl, GEOS_GEOMETRYCOLLECTION, geoms, g2 == NULL ? 1 : 2);
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_reverse(PG_FUNCTION_ARGS)
 */
GEOSGeom ST_Reverse(const GEOSGeom geom)
{
    LWGEOM* lwgeom = GEOS2LWGEOM(geom, 0);
    lwgeom_reverse(lwgeom);
    GEOSGeom ret = LWGEOM2GEOS(lwgeom);
    lwgeom_free(lwgeom);
    return ret;
}

GEOSGeom ST_EndPoint(const GEOSGeometry* geom)
{
    assert (geom && GEOSGeomTypeId_r(hdl, geom) == GEOS_LINESTRING);
    return GEOSGeomGetEndPoint_r(hdl, geom);
}

GEOSGeom ST_Envelope(const GEOSGeom geom)
{
    return GEOSEnvelope_r(hdl, geom);
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_force_clockwise_poly(PG_FUNCTION_ARGS)
 */
GEOSGeom ST_ForceRHR(const GEOSGeom geom)
{
    LWGEOM *lwgeom = GEOS2LWGEOM(geom, 0);

    lwgeom_force_clockwise(lwgeom);

    GEOSGeom ret = LWGEOM2GEOS(lwgeom);

    lwgeom_free(lwgeom);

    return ret;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_makeline(PG_FUNCTION_ARGS)
 */
GEOSGeom ST_MakeLine(const GEOSGeom g1, const GEOSGeom g2)
{
    LWGEOM* lwgeoms[2];
    lwgeoms[0] = GEOS2LWGEOM(g1, 0);
    lwgeoms[1] = GEOS2LWGEOM(g2, 0);

    LWLINE* outline = lwline_from_lwgeom_array(lwgeoms[0]->srid, 2, lwgeoms);

    GEOSGeom ret = LWGEOM2GEOS((LWGEOM *)outline);

    lwgeom_free(lwgeoms[0]);
    lwgeom_free(lwgeoms[1]);

    return ret;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_setpoint_linestring(PG_FUNCTION_ARGS)
 */
GEOSGeom ST_SetPoint(const GEOSGeometry* line, int index, const GEOSGeometry* point)
{
    assert (line && GEOSGeomTypeId(line) == GEOS_LINESTRING);
    assert (index >= 0 && index < GEOSGeomGetNumPoints_r(hdl, line));
    assert (point && GEOSGeomTypeId(point) == GEOS_POINT);

    LWGEOM* lwgeom_line  = GEOS2LWGEOM(line, 0);
    LWGEOM* lwgeom_point = GEOS2LWGEOM(point, 0);

    LWLINE* lwline   = lwgeom_as_lwline(lwgeom_line);       // Both are just casts, no need to free them
    LWPOINT* lwpoint = lwgeom_as_lwpoint(lwgeom_point);     //

    POINT4D newpoint;
    getPoint4d_p(lwpoint->point, 0, &newpoint);

    lwline_setPoint4d(lwline, index, &newpoint);

    GEOSGeom ret = LWGEOM2GEOS(lwgeom_line);

    lwgeom_free(lwgeom_line);
    lwgeom_free(lwgeom_point);

    return ret;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_geos.c):
 *   Datum ST_BuildArea(PG_FUNCTION_ARGS)
 */
GEOSGeom ST_BuildArea(const GEOSGeom geom)
{
    LWGEOM* lwgeom_in = GEOS2LWGEOM(geom, 0);
    LWGEOM* lwgeom_out;

    lwgeom_out = lwgeom_buildarea(lwgeom_in);

    if (!lwgeom_out) {
        return NULL;
    }

    GEOSGeom ret = LWGEOM2GEOS(lwgeom_out);

    lwgeom_free(lwgeom_in) ;
    lwgeom_free(lwgeom_out) ;

    return ret;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_ogc.c):
 *   Datum LWGEOM_geometryn_collection(PG_FUNCTION_ARGS)
 */
GEOSGeom ST_GeometryN(const GEOSGeom geom, int index)
{
    if (index < 0) {
        return NULL;
    }

    LWGEOM* lwgeom = GEOS2LWGEOM(geom, 0);

    LWCOLLECTION* coll = lwgeom_as_lwcollection(lwgeom);
    if (index >= coll->ngeoms) {
        return NULL;
    }

    LWGEOM* subgeom = coll->geoms[index];
    subgeom->srid = coll->srid;

#if 0
    /* COMPUTE_BBOX==TAINTING */
    if ( coll->bbox ) lwgeom_add_bbox(subgeom);
#endif

    GEOSGeom ret = LWGEOM2GEOS(subgeom);

    lwgeom_free(lwgeom);
    lwcollection_free(coll);

    return ret;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_geos_clean.c):
 *   Datum ST_MakeValid(PG_FUNCTION_ARGS)
 */
GEOSGeom ST_MakeValid(const GEOSGeometry* geom)
{
    LWGEOM* lwgeom_in = GEOS2LWGEOM(geom, 0);
    LWGEOM* lwgeom_out = lwgeom_make_valid(lwgeom_in);

    if (!lwgeom_out) {
        lwgeom_free(lwgeom_in);
        return NULL;
    }

    GEOSGeom ret = LWGEOM2GEOS(lwgeom_out);

    lwgeom_free(lwgeom_in);
    lwgeom_free(lwgeom_out);

    return ret;
}

GEOSGeom ST_StartPoint(const GEOSGeometry* geom)
{
    assert (geom && GEOSGeomTypeId_r(hdl, geom) == GEOS_LINESTRING);
    return GEOSGeomGetStartPoint_r(hdl, geom);
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_makepoly(PG_FUNCTION_ARGS)
 */
GEOSGeom ST_MakePolygon(const GEOSGeom geom)
{
    const LWLINE* shell = lwgeom_as_lwline(GEOS2LWGEOM(geom, 0));
    LWPOLY* outpoly = lwpoly_from_lwlines(shell, 0, NULL);

    GEOSGeom ret = LWGEOM2GEOS((LWGEOM *)outpoly);

    lwpoly_free(outpoly);

    return ret;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_closestpoint(PG_FUNCTION_ARGS)
 */
GEOSGeom ST_ClosestPoint(const GEOSGeom g1, const GEOSGeom g2)
{
    LWGEOM* point;
    LWGEOM* lwgeom1 = GEOS2LWGEOM(g1, 0);
    LWGEOM* lwgeom2 = GEOS2LWGEOM(g2, 0);

    point = lw_dist2d_distancepoint(lwgeom1, lwgeom2, lwgeom1->srid, DIST_MIN);

    if (lwgeom_is_empty(point)) {
        return NULL;
    }

    GEOSGeom ret = LWGEOM2GEOS(point);

    lwgeom_free(point);
    lwgeom_free(lwgeom1);
    lwgeom_free(lwgeom2);

    return ret;
}

GEOSGeom ST_CollectionExtract(const GEOSGeometry* geom, int type)
{
    assert (type >= 1 && type <= 3);
    assert (GEOSGeomTypeId(geom) == GEOS_GEOMETRYCOLLECTION);

    vector<GEOSGeom> geoms;
    for (int i = 0; i < GEOSGetNumGeometries_r(hdl, geom); ++i) {
        geoms.push_back(GEOSGeom_clone_r(hdl, GEOSGetGeometryN_r(hdl, geom, i)));
    }

    return GEOSGeom_createCollection_r(
        hdl,
        GEOS_GEOMETRYCOLLECTION,
        geoms.data(),
        geoms.size()
    );
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum ST_RemoveRepeatedPoints(PG_FUNCTION_ARGS)
 */
GEOSGeom ST_RemoveRepeatedPoints(const GEOSGeom geom)
{
    LWGEOM* lwgeom = GEOS2LWGEOM(geom, 0);
    LWGEOM* outgeom = lwgeom_remove_repeated_points(lwgeom);

    GEOSGeom ret = LWGEOM2GEOS(outgeom);

    lwgeom_free(lwgeom);
    lwgeom_free(outgeom);

    return ret;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_ogc.c):
 *   Datum LWGEOM_numpoints_linestring(PG_FUNCTION_ARGS)
 */
int ST_NPoints(const GEOSGeometry* geom)
{
    LWGEOM* lwgeom = GEOS2LWGEOM(geom, 0);

    int ret = lwgeom_count_vertices(lwgeom);

    lwgeom_free(lwgeom);

    return ret;
}

bool bounding_box(const GEOSGeom geom, vector<double>& bbox)
{
    assert (bbox.size() == 0);

    const GEOSCoordSequence* seq = GEOSGeom_getCoordSeq(geom);
    unsigned int size;
    GEOSCoordSeq_getSize(seq, &size);

    if (size < 1) {
        return false;
    }

    double x, y;
    GEOSCoordSeq_getX(seq, 0, &x);
    GEOSCoordSeq_getY(seq, 0, &y);

    double minX = x;
    double minY = y;
    double maxX = x;
    double maxY = y;

    for (int idx = 1; idx < size; ++idx) {
        GEOSCoordSeq_getX(seq, idx, &x);
        GEOSCoordSeq_getY(seq, idx, &y);

        minX = min(minX, x);
        minY = min(minY, y);
        maxX = max(maxX, x);
        maxY = max(maxY, y);
    }

    bbox.push_back(minX);
    bbox.push_back(minY);
    bbox.push_back(maxX);
    bbox.push_back(maxY);

    return true;
}

} // namespace cma
