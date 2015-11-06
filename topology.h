#ifndef __CMA_TOPOLOGY_H
#define __CMA_TOPOLOGY_H

#include <limits>
#include <vector>
#include <geos_c.h>
#include <algorithm>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include <st.h>
#include <types.h>
#include <utils.h>

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
typedef boost::geometry::index::rtree< edge_value, boost::geometry::index::rstar<30000> > edge_idx_t;
typedef boost::geometry::index::rtree< node_value, boost::geometry::index::rstar<30000> > node_idx_t;

/**
 * Other useful types.
 */
typedef boost::geometry::model::polygon<point> polygon;
typedef boost::geometry::model::linestring<point> linestring;
typedef boost::geometry::model::multi_linestring<linestring> multi_linestring;

class Topology
{
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

    void add_edge(edge* e);
    void add_node(node* n);
    void add_face(face* f);

    void output_nodes() const;
    void output_edges() const;
    void output_faces() const;

    const node* closest_and_within_node(const GEOSGeometry* geom, double tolerance);
    const edge* closest_and_within_edge(const GEOSGeometry* geom, double tolerance);

private:
    std::vector<node*> _nodes;
    std::vector<edge*> _edges;
    std::vector<face*> _faces;
    std::vector<relation*> _relations;

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

    template<class T>
    bool _is_in(T hay, const std::vector<T>& stack) const;

    template <class T>
    bool _is_null(T& val) const;

    int _ST_AddFaceSplit(int edgeId, int faceId, bool mbrOnly);
    void GetRingEdges(int edgeId, std::vector<int>& ringEdgeIds, int maxEdges=NULLint);
    void _find_links_to_node(int nodeId, std::vector<edge*>& edges, _span_t& pan, bool span, edge* newEdge, bool isclosed);

    template <class IndexType, class Value>
    void _intersects(IndexType* index, const GEOSGeometry* geom, std::vector<int>& ids);
};

/**
 * Convert a GEOS geometry (GEOS_MULTILINESTRING) to a Boost multi_linestring.
 */
void GEOM2BOOSTMLS(const GEOSGeometry* in, multi_linestring& mls);

/**
 * Convert a GEOS geometry (GEOS_LINESTRING or GEOS_LINEARRING) to a Boost linestring.
 */
void GEOM2BOOSTLS(const GEOSGeometry* in, linestring& ls);

template<class T>
bool Topology::_is_in(T hay, const std::vector<T>& stack) const
{
    return (std::find(stack.begin(), stack.end(), hay) != stack.end());
}

template <class T>
bool Topology::_is_null(T& val) const
{
    return (val == std::numeric_limits<T>::max());
}

/**
 * Delete all elements of a pointer vector and empty it.
 */
template <class T>
void Topology::delete_all(std::vector<T*>& v) {
    for (T* t : v) {
        delete t;
    }
    v.empty();
}

/**
 * Return all items intersecting with geom in the given index.
 */
template <class IndexType, class Value>
void Topology::_intersects(IndexType* index, const GEOSGeometry* geom, std::vector<int>& edgeIds) {
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
        point pt(x, y);
        index->query(boost::geometry::index::intersects(pt), back_inserter(results_s));
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
        index->query(boost::geometry::index::intersects(ls), back_inserter(results_s));

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
