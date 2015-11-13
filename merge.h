#ifndef __CMA_MERGE_H
#define __CMA_MERGE_H

#include <topology.h>

namespace cma {

/**
 * Merge two topologies into one.
 *
 * The first topology will be modified and the second one will be emptied.
 * Note that they MUST be entirely independent topologies.
 */
void merge_topologies(Topology& t1, Topology& t2);

} // namespace cma

#endif // __CMA_MERGE_H
