#ifndef __CMA_ZONES_H
#define __CMA_ZONES_H

#include <vector>
#include <geos_c.h>
#include <ogrsf_frmts.h>

#include <pg.h>
#include <types.h>
#include <utils.h>

namespace cma {

extern GEOSContextHandle_t hdl;

void prepare_zones(PG& db, GEOSHelper& geos, const GEOSGeometry* extent, std::vector<zoneInfo*>& zones, int maxdepth=1);
void write_zones(const std::string& filename, const std::vector<zoneInfo*>& zones, bool overwrite = true);

/**
 * Return a geometry (polygon) which represents the entire world.
 * The resulting geometry must be freed by the caller (GEOSGeom_destroy_r).
 */
GEOSGeometry* world_geom();

/**
 * Return a GEOS geometry reprensenting the provided OGR envelope.
 * The resulting geometry must be freed by the caller (GEOSGeom_destroy_r).
 */
GEOSGeometry* OGREnvelope2GEOSGeom(const OGREnvelope& env);

} // namespace cma

#endif // __CMA_ZONES_H
