#ifndef __CMA_ZONES_H
#define __CMA_ZONES_H

#include <array>
#include <vector>
#include <geos_c.h>
#include <ogrsf_frmts.h>

#include <pg.h>
#include <types.h>
#include <utils.h>

namespace cma {

/**
 * This is to store the spliting operattion so it can be merged
 * back in the correct order.
 */
typedef std::array<int, 4> group_t;
typedef std::pair<int, group_t> depth_group_t;

extern GEOSContextHandle_t hdl;

int prepare_zones(
    GEOSHelper& geos,
    const GEOSGeometry* extent,
    std::vector<zone*>& zones,
    std::vector<depth_group_t>& grouping,
    int maxdepth=1,
    int* nextZoneId=nullptr
);

void write_zones(
    const std::string& filename,
    const std::vector<zone*>& zones,
    bool overwrite = true
);

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
