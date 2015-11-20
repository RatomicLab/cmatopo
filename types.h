#ifndef __CMA_TYPES_H
#define __CMA_TYPES_H

#include <set>
#include <limits>
#include <vector>
#include <memory>

#include <geos_c.h>
#include <ogrsf_frmts.h>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#define NULLint std::numeric_limits<int>::max()
#define NULLdbl std::numeric_limits<double>::max()
#define NULLsizet std::numeric_limits<size_t>::max()

namespace cma {

typedef std::vector<GEOSGeometry *> linesV;
typedef std::pair<OGREnvelope, int> zoneInfo;

class geom_container
{
    friend class boost::serialization::access;

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

  private:
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        if (!geom) {
            ar & 0;
            return;
        }

        GEOSWKBWriter* wkbw = GEOSWKBWriter_create();

        size_t size;
        unsigned char* bin = GEOSWKBWriter_write(wkbw, geom, &size);
        ar & size;
        for (int i = 0; i < size; ++i) {
            ar & bin[i];
        }
        GEOSFree(bin);

        GEOSWKBWriter_destroy(wkbw);
    }

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        size_t size;

        ar & size;
        if (size == 0) {
            geom = nullptr;
            return;
        }

        unsigned char* bin = new unsigned char[size];
        for (int i = 0; i < size; ++i) {
            ar & bin[i];
        }

        GEOSWKBReader* wkbr = GEOSWKBReader_create();
        geom = GEOSWKBReader_read(wkbr, bin, size);
        GEOSWKBReader_destroy(wkbr);

        delete [] bin;
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

class edge : public geom_container
{
    friend class boost::serialization::access;

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

    ~edge() {}

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

  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<geom_container>(*this);
        ar & id;
        ar & start_node;
        ar & end_node;
        ar & next_left_edge;
        ar & next_right_edge;
        ar & abs_next_left_edge;
        ar & abs_next_right_edge;
        ar & left_face;
        ar & right_face;
    }
};

class node : public geom_container
{
    friend class boost::serialization::access;

  public:
    node(): geom_container() {}

    node(const node* other, bool clone_geom = true): node(*other, clone_geom) {}

    node(const node& other, bool clone_geom = true):
        geom_container(other, clone_geom),
        id(other.id),
        containing_face(other.containing_face) {};

    ~node() {}
    bool intersects(const GEOSGeometry* geom);

    int id = NULLint;
    int containing_face = NULLint;

  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<geom_container>(*this);
        ar & id;
        ar & containing_face;
    }
};

class face : public geom_container
{
    friend class boost::serialization::access;

  public:
    face(): geom_container() {}

    face(const face* other, bool clone_geom = true): face(*other, clone_geom) {}

    face(const face& other, bool clone_geom = true):
        geom_container(other, clone_geom),
        id(other.id) {};

    ~face() {}

    int id = NULLint;

  private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<geom_container>(*this);
        ar & id;
    }
};

class relation
{ 
    friend class boost::serialization::access;

public:
    int id = NULLint;     // convenience, not in original table
    int topogeo_id = NULLint;
    int layer_id = NULLint;
    int element_id = NULLint;
    int element_type = NULLint;

private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & id;
        ar & topogeo_id;
        ar & layer_id;
        ar & element_id;
        ar & element_type;
    }
};

typedef std::set<int> edgeid_set;
typedef std::shared_ptr<edgeid_set> edgeid_set_ptr;

} // namespace cma

#endif // __CMA_TYPES_H
