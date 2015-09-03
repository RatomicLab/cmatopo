#ifndef __CMA_ST_H
#define __CMA_ST_H

#include <types.h>

#include <limits>
#include <string>
#include <vector>

namespace cma {

/**
 * Port of some PostGIS/topology functions.
 */

bool ST_Equals(const GEOSGeom g1, const GEOSGeom g2);
bool ST_DWithin(const GEOSGeom g1, const GEOSGeom g2, double tolerance);
bool ST_IsEmpty(GEOSContextHandle_t geos, const GEOSGeom geom);
bool ST_Contains(const GEOSGeom g1, const GEOSGeom g2);
bool ST_OrderingEquals(const GEOSGeom g1, const GEOSGeom g2);

double ST_X(const GEOSGeom geom);
double ST_Y(const GEOSGeom geom);
double ST_Distance(const GEOSGeom g1, const GEOSGeom g2);
double _ST_MinTolerance(const GEOSGeom geom);
double ST_Azimuth(GEOSContextHandle_t geos, const GEOSGeometry* g1, const GEOSGeometry* g2);

GEOSGeom ST_Snap(const GEOSGeom g1, const GEOSGeom g2, double tolerance);
GEOSGeom ST_Split(const GEOSGeometry* in, const GEOSGeometry* blade_in);
GEOSGeom ST_PointN(GEOSContextHandle_t geos, const GEOSGeometry* line, int index);
GEOSGeom ST_Collect(GEOSGeometry* g1, GEOSGeometry* g2 = NULL);
GEOSGeom ST_Reverse(const GEOSGeom geom);
GEOSGeom ST_EndPoint(GEOSContextHandle_t geos, const GEOSGeometry* geom);
GEOSGeom ST_Envelope(const GEOSGeom geom);
GEOSGeom ST_ForceRHR(const GEOSGeom geom);
GEOSGeom ST_MakeLine(const GEOSGeom g1, const GEOSGeom g2);
GEOSGeom ST_SetPoint(GEOSContextHandle_t _geos, const GEOSGeometry* line, int index, const GEOSGeometry* point);
GEOSGeom ST_BuildArea(const GEOSGeom geom);
GEOSGeom ST_GeometryN(const GEOSGeom geom, int index);
GEOSGeom ST_MakeValid(const GEOSGeometry* geom);
GEOSGeom ST_StartPoint(GEOSContextHandle_t geos, const GEOSGeometry* geom);
GEOSGeom ST_MakePolygon(const GEOSGeom geom);
GEOSGeom ST_ClosestPoint(const GEOSGeom g1, const GEOSGeom g2);
GEOSGeom ST_CollectionExtract(GEOSContextHandle_t _geos, const GEOSGeometry* collection, int type);
GEOSGeom ST_RemoveRepeatedPoints(const GEOSGeom geom);

int ST_NPoints(const GEOSGeometry* geom);

bool bounding_box(const GEOSGeom geom, std::vector<double>& bbox);

/**
 * Find, if it exists, the geometry from a set of geometries (others) which is the closest within
 * a specified tolerance.
 */
template <class T>
const T* closest_and_within(const GEOSGeom geom, const std::vector<T*>& others, double tolerance)
{
    const T* item = NULL;
    double previousDistance = std::numeric_limits<double>::max();
    for (const T* other : others) {
        if (ST_DWithin(other->geom, geom, tolerance) && (previousDistance = ST_Distance(geom, other->geom)) < previousDistance) {
            item = other;
        }
    }
    return item;
}

} // namespace cma

#endif // __CMA_ST_H
