#ifndef __CMA_MERGE_H
#define __CMA_MERGE_H

#include <topology.h>

#include <vector>

#include <zones.h>

namespace cma {

class zone;

typedef enum {
    ABOVE,
    BELOW,
    RIGHT,
    LEFT,
    OTHER
} direction_type;

/**
 * Merge two topologies into one.
 *
 * The first topology will be modified and the second one will be emptied.
 * Note that they MUST be entirely independent topologies.
 */
void merge_topologies(Topology& t1, Topology& t2);

/**
 * Proceed with a pair-wise merge and return merged topologies
 * as well as merged zones.
 */
void merge_topologies(
    const std::vector<zone*>& zones,
    std::vector<Topology*>& topologies,
    std::vector<zone*>& new_zones,
    std::vector<Topology*>& new_topologies
);

/**
 * Get the next 4-grouped zones that can be merged.
 */
void get_next_groups(
    std::vector<depth_group_t> all_groups,
    std::vector<depth_group_t> next_groups
);

/**
 * Get the width of an envelope.
 */
double width(const OGREnvelope& envelope);

/**
 * Get the height of an envelope.
 */
double height(const OGREnvelope& envelope);

direction_type position(const OGREnvelope& e1, const OGREnvelope& e2);

} // namespace cma

#endif // __CMA_MERGE_H
