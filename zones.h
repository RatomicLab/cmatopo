#ifndef __CMA_ZONES_H
#define __CMA_ZONES_H

#include <vector>
#include <ogrsf_frmts.h>

#include <types.h>

namespace cma {

void prepare_zones(const linesV& lines, const OGREnvelope& extent, std::vector<zoneInfo*>& zones, int maxdepth=1);
void write_zones(const std::string& filename, const std::vector<zoneInfo*>& zones, bool overwrite = true);

} // namespace cma

#endif // __CMA_ZONES_H
