#ifndef __CMA_TOPOLOGY_H
#define __CMA_TOPOLOGY_H

#include <set>
#include <limits>
#include <memory>
#include <vector>
#include <geos_c.h>
#include <algorithm>

#include <boost/tuple/tuple.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include <st.h>
#include <types.h>
#include <utils.h>
#include <transaction.h>

#define DEFAULT_TOLERANCE 1.0

namespace cma {

/**
 * Internal structure used for _find_links_to_node azimuth computation.
 */
class _span_t {
  public:
    double myaz = NULLdbl;
    double minaz = NULLdbl;
    double maxaz = NULLdbl;
    int nextCW = NULLint;
    int nextCCW = NULLint;
    bool was_isolated = false;
};

/**
 * Types used for indexing.
 */
typedef boost::geometry::model::d2::point_xy<double> point;
typedef boost::geometry::model::box<point> box;
typedef std::pair<box,   int> edge_value;
typedef std::pair<point, int> node_value;
typedef boost::geometry::index::rtree< edge_value, boost::geometry::index::rstar<100000> > edge_idx_t;
typedef boost::geometry::index::rtree< node_value, boost::geometry::index::rstar<100000> > node_idx_t;

/**
 * Other useful types.
 */
typedef boost::geometry::model::polygon<point> polygon;
typedef boost::geometry::model::linestring<point> linestring;
typedef boost::geometry::model::multi_linestring<linestring> multi_linestring;

class Topology
{
    friend class FaceTransaction;
    friend class EdgeTransaction;
    friend class NodeTransaction;
    friend class TopologyTransaction;
    friend class AddFaceIndexTransaction;
    friend class RemoveFaceIndexTransaction;
    friend class AddRelationTransaction;

    friend void merge_topologies(Topology&, Topology&);

public:
    Topology(GEOSHelper& geos);
    ~Topology();

    template <class T>
    void delete_all(std::vector<T*>& v);

    /**
     * Add an edge (and it's endpoints) to the topology.
     */
    void TopoGeo_AddLineString(GEOSGeom line, std::vector<int>& edgeIds, double tolerance=0.);
    int ST_AddEdgeModFace(int start_node, int end_node, GEOSGeometry* geom);

    /*****************/

    int TopoGeo_AddPoint(GEOSGeom point, double tolerance=0.);

    const node* get_node_at(double x, double y) {
        std::vector<node_value> nodes;
        _node_idx->query(boost::geometry::index::intersects(point(x,y)), back_inserter(nodes));
        assert (nodes.size() == 0 || nodes.size() == 1);

        if (nodes.size() == 0) {
            return NULL;
        }

        return _nodes[nodes[0].second];
    }

    int ST_ChangeEdgeGeom(int edgeId, const GEOSGeom point);
    int ST_ModEdgeSplit(int edgeId, const GEOSGeom point);
    void _ST_AdjacentEdges(int nodeId, int edgeId, std::vector<int>& edges);
    void GetNodeEdges(int nodeId, std::vector<int>& edgeIds);
    GEOSGeom ST_GetFaceGeometry(int faceId);
    int ST_AddIsoNode(int faceId, const GEOSGeom point);

    void commit();
    void rollback();

    void output() const;
    void output_nodes() const;
    void output_edges() const;
    void output_faces() const;
    void output_relations() const;

    const node* closest_and_within_node(const GEOSGeometry* geom, double tolerance);
    const edge* closest_and_within_edge(const GEOSGeometry* geom, double tolerance);

private:
    std::vector<node*> _nodes;
    std::vector<edge*> _edges;
    std::vector<face*> _faces;
    std::vector<relation*> _relations;

    /**
     * Rollback data members
     */
    std::vector<TopologyTransaction*>* _transactions = nullptr;
    std::vector<int>* _inserted_nodes = nullptr;
    std::vector<int>* _inserted_edges = nullptr;
    std::vector<int>* _inserted_faces = nullptr;

    /**
     * Temporary vector to track geometry deletion
     * when a commit/rollback operation occurs.
     */
    std::set<GEOSGeometry*>* _tr_track_geom = nullptr;

    /**
     * GEOS helper class.
     */
    GEOSHelper& _geos;

    /**
     * Edge envelope index.
     */
    edge_idx_t* _edge_idx = NULL;

    /**
     * Edge envelope + tolerance around it index.
     */
    edge_idx_t* _edge_tol_idx = NULL;

    /**
     * Node index.
     */
    node_idx_t* _node_idx = NULL;

    /**
     * Node + tolerance around it.
     * edge_idx_t* is not a typo, we store boxes.
     */
    edge_idx_t* _node_tol_idx = NULL;

    /**
     * Total linestrings that were added to this topology.
     */
    uint64_t _totalCount = 0;

    /**
     * Index of edges left faces.
     */
    std::vector<edgeid_set_ptr>* _left_faces_idx = nullptr;

    /**
     * Index of edges right faces
     */
    std::vector<edgeid_set_ptr>* _right_faces_idx = nullptr;

    /**
     * Face geometry cache.
     */
    std::vector<GEOSGeometry*>* _face_geometries = nullptr;

    /**
     * Temporary vector for ST_GetFaceGeometry operations.
     */
    std::vector<GEOSGeometry*>* _gfg_geometries = nullptr;

    void add_edge(edge* e);
    void add_node(node* n);
    void add_face(face* f);

    void remove_edge(int edgeId);
    void remove_node(int nodeId);
    void remove_face(int faceId);

    void _update_left_face(edge* e, int faceId);
    void _update_right_face(edge* e, int faceId);

    void _face_edges(int faceId, edgeid_set& edges);

    void _empty(bool free_items=true);

    int _ST_AddFaceSplit(int edgeId, int faceId, bool mbrOnly);
    void GetRingEdges(int edgeId, std::vector<int>& ringEdgeIds, int maxEdges=NULLint);
    void _find_links_to_node(int nodeId, std::vector<edge*>& edges, _span_t& pan, bool span, edge* newEdge, bool isclosed);

    template <class IndexType, class Value>
    void _intersects(IndexType* index, const GEOSGeometry* geom, std::vector<int>& ids, double tolerance = 0.);
};

/**
 * Convert a GEOS geometry (GEOS_MULTILINESTRING) to a Boost multi_linestring.
 */
void GEOM2BOOSTMLS(const GEOSGeometry* in, multi_linestring& mls);

/**
 * Convert a GEOS geometry (GEOS_LINESTRING or GEOS_LINEARRING) to a Boost linestring.
 */
void GEOM2BOOSTLS(const GEOSGeometry* in, linestring& ls);

/**
 * Delete all elements of a pointer vector and empty it.
 */
template <class T>
void Topology::delete_all(std::vector<T*>& v) {
    for (T* t : v) {
        if (!t) continue;
        delete t;
    }
    v.clear();
}

/**
 * Return all items intersecting with geom in the given index.
 */
template <class IndexType, class Value>
void Topology::_intersects(IndexType* index, const GEOSGeometry* geom, std::vector<int>& edgeIds, double tolerance) {
    assert (index);
    assert (geom);
    assert (edgeIds.size() == 0);

    std::vector<Value> results_s;

    switch (GEOSGeomTypeId_r(hdl, geom))
    {
    case GEOS_POINT: {
        double x, y;
        GEOSGeomGetX_r(hdl, geom, &x);
        GEOSGeomGetY_r(hdl, geom, &y);

        if (tolerance == 0.) {
            point pt(x, y);
            index->query(boost::geometry::index::intersects(pt), back_inserter(results_s));
        }
        else {
            box env;

            env.min_corner().x(x - tolerance/2);
            env.min_corner().y(y - tolerance/2);
            env.max_corner().x(x + tolerance/2);
            env.max_corner().y(x + tolerance/2);

            index->query(boost::geometry::index::intersects(env), back_inserter(results_s));
        }

        break;
    }
    case GEOS_LINESTRING:
    case GEOS_POLYGON: {
        const GEOSGeometry* g = geom;
        if (GEOSGeomTypeId_r(hdl, geom) == GEOS_POLYGON) {
            g = GEOSGetExteriorRing_r(hdl, geom);
        }

        linestring ls;
        GEOM2BOOSTLS(g, ls);

        box env;
        boost::geometry::envelope(ls, env);

        if (tolerance > 0.) {
            env.min_corner().x(env.min_corner().x() - tolerance/2);
            env.min_corner().y(env.min_corner().y() - tolerance/2);
            env.max_corner().x(env.max_corner().x() + tolerance/2);
            env.max_corner().y(env.max_corner().y() + tolerance/2);
        }

        index->query(boost::geometry::index::intersects(env), back_inserter(results_s));

        break;
    }
    case GEOS_MULTILINESTRING: {
        multi_linestring mls;
        GEOM2BOOSTMLS(geom, mls);
        index->query(boost::geometry::index::intersects(mls), back_inserter(results_s));

        break;
    }
    default:
        std::cout << GEOSGeomTypeId_r(hdl, geom)  << endl;
        throw std::invalid_argument("unsupported geometry type.");
    };

    // Sort by id. Strickly speaking, this is not necessary but we want to
    // maintain order to compare our results with PostGIS.
    sort(results_s.begin(), results_s.end(), [](const Value& a, const Value& b) {
        return a.second < b.second;
    });

    if (results_s.size() > 0) {
        edgeIds.reserve(results_s.size()+1);
        transform(results_s.begin(), results_s.end(), back_inserter(edgeIds), [](const Value& a) {
            return a.second;
        });
    }
}

} // namespace cma

#endif // __CMA_TOPOLOGY_H
