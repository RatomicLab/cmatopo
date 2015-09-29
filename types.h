#ifndef __CMA_TYPES_H
#define __CMA_TYPES_H

#include <limits>
#include <vector>

#include <geos_c.h>
#include <ogrsf_frmts.h>

#define NULLint std::numeric_limits<int>::max()
#define NULLdbl std::numeric_limits<double>::max()
#define NULLsizet std::numeric_limits<size_t>::max()

namespace cma {

typedef std::vector<GEOSGeometry *> linesV;
typedef std::pair<OGREnvelope, int> zoneInfo;

class edge {
  public:
    edge() {}

    edge(const edge* other): edge(*other) {}

    edge(const edge& other):
        id(other.id),
        start_node(other.start_node),
        end_node(other.end_node),
        next_left_edge(other.next_left_edge),
        abs_next_left_edge(other.abs_next_left_edge),
        next_right_edge(other.next_right_edge),
        abs_next_right_edge(other.abs_next_right_edge),
        left_face(other.left_face),
        right_face(other.right_face),
        prev_left_edge(other.prev_left_edge),
        prev_right_edge(other.prev_right_edge),
        geom(GEOSGeom_clone(other.geom)) {}

    ~edge();

    int id = NULLint;
    int start_node = NULLint;
    int end_node   = NULLint;
    int next_left_edge  = NULLint;
    int next_right_edge = NULLint;
    int abs_next_left_edge  = NULLint;
    int abs_next_right_edge = NULLint;
    int left_face  = NULLint;
    int right_face = NULLint;
    int prev_left_edge  = NULLint;       // convenience
    int prev_right_edge = NULLint;      // convenience
    GEOSGeom geom = NULL;

    bool intersects(const GEOSGeometry* geom);

  private:
    GEOSGeometry* _envelope = NULL;
    const GEOSGeometry* envelope();
};

class node {
  public:
    ~node();

    int id = NULLint;
    int containing_face = NULLint;
    GEOSGeom geom = NULL;

    bool intersects(const GEOSGeometry* geom);

  private:
    GEOSGeometry* _envelope = NULL;
    const GEOSGeometry* envelope();
};

class face {
  public:
    ~face();

    int id = NULLint;
    GEOSGeom mbr = NULL;

    bool intersects(const GEOSGeometry* geom);

  private:
    GEOSGeometry* _envelope = NULL;
    const GEOSGeometry* envelope();
};

typedef struct {
    int element_id = NULLint;
    int element_type = NULLint;
} relation;

inline void sort_zones_by_line_count(std::vector<zoneInfo*>& zones, bool reverse_sort=false) {
    std::sort(zones.begin(), zones.end(), [](zoneInfo* a, zoneInfo* b) {
        return a < b;
    });
    if (reverse_sort) {
        reverse(zones.begin(), zones.end());
    }
}

} // namespace cma

#endif // __CMA_TYPES_H
