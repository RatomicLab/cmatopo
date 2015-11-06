#ifndef __CMA_ST_H
#define __CMA_ST_H

#include <types.h>

#include <limits>
#include <string>
#include <vector>

namespace cma {

extern GEOSContextHandle_t hdl;

/**
 * Port of some PostGIS/topology functions.
 */

bool ST_Equals(const GEOSGeom g1, const GEOSGeom g2);
bool ST_DWithin(const GEOSGeometry* g1, const GEOSGeometry* g2, double tolerance);
bool ST_IsEmpty(const GEOSGeom geom);
bool ST_Contains(const GEOSGeom g1, const GEOSGeom g2);
bool ST_OrderingEquals(const GEOSGeom g1, const GEOSGeom g2);

double ST_X(const GEOSGeom geom);
double ST_Y(const GEOSGeom geom);
double ST_Azimuth(const GEOSGeometry* g1, const GEOSGeometry* g2);
double ST_Distance(const GEOSGeometry* g1, const GEOSGeometry* g2);
double _ST_MinTolerance(const GEOSGeom geom);

GEOSGeom ST_Snap(const GEOSGeom g1, const GEOSGeom g2, double tolerance);
GEOSGeom ST_Split(const GEOSGeometry* in, const GEOSGeometry* blade_in);
GEOSGeom ST_PointN(const GEOSGeometry* line, int index);
GEOSGeom ST_Collect(GEOSGeometry* g1, GEOSGeometry* g2 = NULL);
GEOSGeom ST_Reverse(const GEOSGeom geom);
GEOSGeom ST_AddPoint(GEOSGeometry* line, GEOSGeometry* pt, int where = -1);
GEOSGeom ST_EndPoint(const GEOSGeometry* geom);
GEOSGeom ST_Envelope(const GEOSGeom geom);
GEOSGeom ST_ForceRHR(const GEOSGeom geom);
GEOSGeom ST_MakeLine(const GEOSGeom g1, const GEOSGeom g2);
GEOSGeometry* ST_MakeLine(const std::vector<const GEOSGeometry*>& geometries);
GEOSGeom ST_SetPoint(const GEOSGeometry* line, int index, const GEOSGeometry* point);
GEOSGeom ST_BuildArea(const GEOSGeom geom);
GEOSGeom ST_GeometryN(const GEOSGeom geom, int index);
GEOSGeom ST_MakeValid(const GEOSGeometry* geom);
GEOSGeom ST_StartPoint(const GEOSGeometry* geom);
GEOSGeom ST_MakePolygon(const GEOSGeom geom);
GEOSGeom ST_ClosestPoint(const GEOSGeom g1, const GEOSGeom g2);
GEOSGeom ST_CollectionExtract(GEOSGeometry* collection, int type);
GEOSGeom ST_RemoveRepeatedPoints(const GEOSGeom geom);

int ST_NPoints(const GEOSGeometry* geom);

bool bounding_box(const GEOSGeom geom, std::vector<double>& bbox);
bool is_collection(const GEOSGeometry* geom);

} // namespace cma

#endif // __CMA_ST_H
