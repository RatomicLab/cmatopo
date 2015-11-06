#include <st.h>

#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>

extern "C" {
    #include <liblwgeom.h>
    #include <lwgeom_geos.h>
}

using namespace std;

namespace cma {

const float DIST_MIN = 1.;

bool ST_Equals(const GEOSGeom g1, const GEOSGeom g2)
{
    return GEOSWithin_r(hdl, g1, g2) == 1 && GEOSWithin_r(hdl, g2, g1) == 1;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_dwithin(PG_FUNCTION_ARGS)
 */
bool ST_DWithin(const GEOSGeometry* g1, const GEOSGeometry* g2, double tolerance)
{
    assert (tolerance >= 0.);

    double distance;

    GEOSDistance_r(hdl, g1, g2, &distance);

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
    return GEOSContains_r(hdl, g1, g2);
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_same(PG_FUNCTION_ARGS)
 */
bool ST_OrderingEquals(const GEOSGeom g1, const GEOSGeom g2)
{
    if (GEOSGeomTypeId_r(hdl, g1) != GEOSGeomTypeId_r(hdl, g2)) {
        return false;
    }

    return GEOSEqualsExact_r(hdl, g1, g2, 0.) == 1;
}

double ST_X(const GEOSGeom geom)
{
    double x;
    GEOSGeomGetX_r(hdl, geom, &x);
    return x;
}

double ST_Y(const GEOSGeom geom)
{
    double y;
    GEOSGeomGetY_r(hdl, geom, &y);
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
double ST_Distance(const GEOSGeometry* g1, const GEOSGeometry* g2)
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
    vector<double> fbbox;
    transform(bbox.begin(), bbox.end(), back_inserter(fbbox), [](double c) {
        return fabs(c);
    });

    double _max = *max_element(fbbox.begin(), fbbox.end());
    double val = _max != 0. ? _max : 1.;

    return 3.6 * pow(10,  -(15 - log10(val)));
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

    GEOSGeom ret = LWGEOM2GEOS(lwgeom_out);

    lwgeom_free(lwblade_in);
    lwgeom_free(lwgeom_in);
    lwgeom_free(lwgeom_out);

    return ret;
}

GEOSGeom ST_PointN(const GEOSGeometry* line, int index)
{
    assert (line && GEOSGeomTypeId_r(hdl, line) == GEOS_LINESTRING);
    assert (index > 0 && index <= GEOSGeomGetNumPoints_r(hdl, line)); // index is one-based in PostGIS
    return GEOSGeomGetPointN_r(hdl, line, index-1);
}

GEOSGeom ST_Collect(GEOSGeometry* g1, GEOSGeometry* g2)
{
    GEOSGeometry* geoms[] = {g1, g2};

    // After this point, g1 and g2 should not be freed.
    // Only the created collection should.

    return GEOSGeom_createCollection_r(hdl, GEOS_MULTILINESTRING, geoms, g2 == NULL ? 1 : 2);
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

GEOSGeom ST_AddPoint(GEOSGeometry* line, GEOSGeometry* point, int where)
{
    assert (line);
    assert (point);
    assert (GEOSGeomTypeId_r(hdl, line) == GEOS_LINESTRING);
    assert (GEOSGeomTypeId_r(hdl, point) == GEOS_POINT);
    assert (where >= -1);

    LWGEOM* lwgeom_line = GEOS2LWGEOM(line, 0);
    LWGEOM* lwgeom_point = GEOS2LWGEOM(point, 0);

    LWLINE* lwline = lwgeom_as_lwline(lwgeom_line);
    LWPOINT* lwpoint = lwgeom_as_lwpoint(lwgeom_point);

    if (where == -1) {
        where = lwline->points->npoints;
    }

    LWLINE* linecopy = lwgeom_as_lwline(lwgeom_clone_deep(lwgeom_line));
    lwline_add_lwpoint(linecopy, lwpoint, where);

    GEOSGeometry* ret = LWGEOM2GEOS(lwline_as_lwgeom(linecopy));

    lwline_free(linecopy);
    lwgeom_free(lwgeom_line);
    lwgeom_free(lwgeom_point);

    return ret;
}

GEOSGeom ST_EndPoint(const GEOSGeometry* geom)
{
    assert (geom && GEOSGeomTypeId_r(hdl, geom) == GEOS_LINESTRING);
    return GEOSGeomGetEndPoint_r(hdl, geom);
}

GEOSGeom ST_Envelope(const GEOSGeom geom)
{
    LWGEOM *lwgeom = GEOS2LWGEOM(geom, 0);
    GBOX box;
    lwgeom_calculate_gbox(lwgeom, &box);

    if ((box.xmin == box.xmax) && (box.ymin == box.ymax) ||
        (box.xmin == box.xmax) || (box.ymin == box.ymax))
    {
        lwgeom_free(lwgeom);
        return GEOSEnvelope_r(hdl, geom);
    }

    // taken from PostGIS: lwgeom_functions_basic.c
    // we don't use GEOSEnvelope_r to get the same point
    // order than PostGIS.

    POINT4D pt;
    LWPOLY *poly;
    POINTARRAY *pa;

    POINTARRAY **ppa = (POINTARRAY **)lwalloc(sizeof(POINTARRAY*));
    pa = ptarray_construct_empty(0, 0, 5);
    ppa[0] = pa;

    /* Assign coordinates to POINT2D array */
    pt.x = box.xmin;
    pt.y = box.ymin;
    ptarray_append_point(pa, &pt, LW_TRUE);
    pt.x = box.xmin;
    pt.y = box.ymax;
    ptarray_append_point(pa, &pt, LW_TRUE);
    pt.x = box.xmax;
    pt.y = box.ymax;
    ptarray_append_point(pa, &pt, LW_TRUE);
    pt.x = box.xmax;
    pt.y = box.ymin;
    ptarray_append_point(pa, &pt, LW_TRUE);
    pt.x = box.xmin;
    pt.y = box.ymin;
    ptarray_append_point(pa, &pt, LW_TRUE);

    /* Construct polygon  */
    poly = lwpoly_construct(GEOSGetSRID_r(hdl, geom), NULL, 1, ppa);
    GEOSGeometry* ret = LWGEOM2GEOS(lwpoly_as_lwgeom(poly));
    lwpoly_free(poly);
    lwgeom_free(lwgeom);

    return ret;
    // return GEOSEnvelope_r(hdl, geom);
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
    vector<const GEOSGeometry*> geometries;
    geometries.push_back(g1);
    geometries.push_back(g2);

    return ST_MakeLine(geometries);
}

GEOSGeometry* ST_MakeLine(const std::vector<const GEOSGeometry*>& geometries)
{
    assert (geometries.size() > 0);

    LWGEOM** lwgeoms = new LWGEOM*[geometries.size()];

    int i = 0;
    for (const GEOSGeometry* g : geometries) {
        int t = GEOSGeomTypeId_r(hdl ,g);
        assert (t == GEOS_POINT || t == GEOS_LINESTRING || t == GEOS_LINEARRING);
        lwgeoms[i++] = GEOS2LWGEOM(g, 0);
    }

    LWLINE* outline = lwline_from_lwgeom_array(lwgeoms[0]->srid, geometries.size(), lwgeoms);

    GEOSGeometry* ret = LWGEOM2GEOS((LWGEOM *)outline);

    for (i = 0; i < geometries.size(); ++i) {
        lwgeom_free(lwgeoms[i]);
    }

    delete [] lwgeoms;

    return ret;
}

/**
 * Mostly equivalent to the following function (postgis/lwgeom_functions_basic.c):
 *   Datum LWGEOM_setpoint_linestring(PG_FUNCTION_ARGS)
 */
GEOSGeom ST_SetPoint(const GEOSGeometry* line, int index, const GEOSGeometry* point)
{
    // index is 0-based in PostGIS
    assert (line && GEOSGeomTypeId_r(hdl, line) == GEOS_LINESTRING);
    assert (index >= 0 && index < GEOSGeomGetNumPoints_r(hdl, line));
    assert (point && GEOSGeomTypeId_r(hdl, point) == GEOS_POINT);

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
        lwgeom_free(lwgeom_in);
        return NULL;
    }

    GEOSGeom ret = LWGEOM2GEOS(lwgeom_out);

    lwgeom_free(lwgeom_in);
    lwgeom_free(lwgeom_out);

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
    return GEOSInterpolate_r(hdl, g1, GEOSProject_r(hdl, g1, g2));
}

/**
 * Extract geometries of the specified type from the collection.
 *
 * Supported types:
 *  - point
 *  - linestring
 *  - linearring
 *  - polygon
 *
 * If the input geometry is not a collection and of the requested type,
 * it is returned as-is (not a copy). If it is not a collection and not of
 * the requested type, an empty geometry of the specified type is returned.
 */
GEOSGeometry* ST_CollectionExtract(GEOSGeometry* geom, int type)
{
    assert (type >= 1 && type <= 4);

    if (!is_collection(geom)) {
        if (GEOSGeomTypeId_r(hdl, geom) == type) {
            return geom;
        }

        GEOSGeometry* ret = NULL;
        switch (type)
        {
            case GEOS_POINT:
                ret = GEOSGeom_createEmptyPoint_r(hdl);
                break;
            case GEOS_LINESTRING:
                ret = GEOSGeom_createEmptyLineString_r(hdl);
                break;
            case GEOS_LINEARRING:
                ret = GEOSGeom_createLinearRing_r(hdl, NULL);   // this has not been tested and might not work
                break;
            case GEOS_POLYGON:
                ret = GEOSGeom_createEmptyPolygon_r(hdl);
                break;
        };
        assert (ret != NULL);
        return ret;
    }

    vector<GEOSGeom> geoms;
    for (int i = 0; i < GEOSGetNumGeometries_r(hdl, geom); ++i) {
        if (GEOSGeomTypeId_r(hdl, GEOSGetGeometryN_r(hdl, geom, i)) == type) {
            geoms.push_back(GEOSGeom_clone_r(hdl, GEOSGetGeometryN_r(hdl, geom, i)));
        }
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

    const GEOSCoordSequence* seq = GEOSGeom_getCoordSeq_r(hdl, geom);
    unsigned int size;
    GEOSCoordSeq_getSize_r(hdl, seq, &size);

    if (size < 1) {
        return false;
    }

    double x, y;
    GEOSCoordSeq_getX_r(hdl, seq, 0, &x);
    GEOSCoordSeq_getY_r(hdl, seq, 0, &y);

    double minX = x;
    double minY = y;
    double maxX = x;
    double maxY = y;

    for (int idx = 1; idx < size; ++idx) {
        GEOSCoordSeq_getX_r(hdl, seq, idx, &x);
        GEOSCoordSeq_getY_r(hdl, seq, idx, &y);

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

bool is_collection(const GEOSGeometry* geom)
{
    int type = GEOSGeomTypeId_r(hdl, geom);
    return (type == GEOS_GEOMETRYCOLLECTION
        || type == GEOS_MULTIPOINT
        || type == GEOS_MULTILINESTRING
        || type == GEOS_MULTIPOLYGON);
}

} // namespace cma
