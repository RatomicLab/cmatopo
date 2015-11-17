#ifndef __CMA_TYPES_H
#define __CMA_TYPES_H

#include <set>
#include <limits>
#include <vector>
#include <memory>

#include <geos_c.h>
#include <ogrsf_frmts.h>

#define NULLint std::numeric_limits<int>::max()
#define NULLdbl std::numeric_limits<double>::max()
#define NULLsizet std::numeric_limits<size_t>::max()

namespace cma {

typedef std::vector<GEOSGeometry *> linesV;
typedef std::pair<OGREnvelope, int> zoneInfo;

class geom_container
{
  public:
      geom_container() {}
      geom_container(const geom_container& other, bool clone = true): geom(other.geom) {
          if (clone) {
              geom = GEOSGeom_clone(other.geom);
          }
      }

      virtual ~geom_container();

      const GEOSGeometry* envelope();
      const GEOSPreparedGeometry* prepared();

      virtual bool intersects(const GEOSGeometry* geom);

      GEOSGeometry* geom = NULL;

  protected:
      GEOSGeometry* _envelope = NULL;
      const GEOSPreparedGeometry* _prepared = NULL;
};

class edge : public geom_container {
  public:
    edge(): geom_container() {}

    edge(const edge* other, bool clone_geom = true): edge(*other, clone_geom) {}

    edge(const edge& other, bool clone_geom = true):
        geom_container(other, clone_geom),
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
        prev_right_edge(other.prev_right_edge) {};

    virtual ~edge() {}

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
    int prev_right_edge = NULLint;       // convenience
};

class node : public geom_container {
  public:
    node(): geom_container() {}

    node(const node* other, bool clone_geom = true): node(*other, clone_geom) {}

    node(const node& other, bool clone_geom = true):
        geom_container(other, clone_geom),
        id(other.id),
        containing_face(other.containing_face) {};

    virtual ~node() {}
    virtual bool intersects(const GEOSGeometry* geom);

    int id = NULLint;
    int containing_face = NULLint;
};

class face : public geom_container {
  public:
    face(): geom_container() {}

    face(const face* other, bool clone_geom = true): face(*other, clone_geom) {}

    face(const face& other, bool clone_geom = true):
        geom_container(other, clone_geom),
        id(other.id) {};

    virtual ~face() {}

    int id = NULLint;
};

typedef struct {
    int id = NULLint;     // convenience, not in original table
    int topogeo_id = NULLint;
    int layer_id = NULLint;
    int element_id = NULLint;
    int element_type = NULLint;
} relation;

typedef std::set<int> edgeid_set;
typedef std::shared_ptr<edgeid_set> edgeid_set_ptr;

} // namespace cma

#endif // __CMA_TYPES_H
