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

class Topology
{
public:
    typedef boost::geometry::model::d2::point_xy<double> point;
    typedef boost::geometry::model::box<point> box;
    //typedef boost::geometry::model::polygon<point> poly;
    typedef std::pair<box,   int> edge_value;
    //typedef std::pair<poly,  int> edge_value;
    typedef std::pair<point, int> node_value;
    typedef boost::geometry::index::rtree< edge_value, boost::geometry::index::rstar<30000> > edge_idx_t;
    typedef boost::geometry::index::rtree< node_value, boost::geometry::index::rstar<30000> > node_idx_t;

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
        for (node* n : _nodes) {
            if (!n) continue;

            // TODO: speedup by testing if n->geom bounding box contains the point POINT(x,y)
            if (ST_X(n->geom) == x && ST_Y(n->geom) == y) {
                return n;
            }
        }
        return NULL;
    }

    int ST_ChangeEdgeGeom(int edgeId, const GEOSGeom point);
    int ST_ModEdgeSplit(int edgeId, const GEOSGeom point);
    void _ST_AdjacentEdges(int nodeId, int edgeId, std::vector<int>& edges);
    void GetNodeEdges(int nodeId, std::vector<int>& edgeIds);
    GEOSGeom ST_GetFaceGeometry(int faceId);
    int ST_AddIsoNode(int faceId, const GEOSGeom point);

    void add_edge(edge* e);
    void add_node(node* n);

    void output_nodes() const;
    void output_edges() const;

    const edge* closest_and_within_edge(const GEOSGeometry* geom, double tolerance);

    template <class T>
    bool contains(const edge* a, const T* b);

    template <class T>
    bool intersects(const edge* a, const T* b);

private:
    std::vector<node*> _nodes;
    std::vector<edge*> _edges;
    std::vector<face*> _faces;
    std::vector<relation*> _relations;

    GEOSHelper& _geos;

    edge_idx_t* _edge_idx = NULL;
    node_idx_t* _node_idx = NULL;

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
 * Detect if a's bounding box contains b's bounding box, using the edge index.
 */
template <class T>
bool Topology::contains(const edge* a, const T* b) {
    const GEOSGeometry* bshell = GEOSGetExteriorRing_r(hdl, b->envelope());
    GEOSGeometry* sp = GEOSGeomGetPointN_r(hdl, bshell, 0);
    GEOSGeometry* ep = GEOSGeomGetPointN_r(hdl, bshell, 3);
    double spX, spY, epX, epY;
    GEOSGeomGetX_r(hdl, sp, &spX);
    GEOSGeomGetY_r(hdl, sp, &spY);
    GEOSGeomGetX_r(hdl, ep, &epX);
    GEOSGeomGetY_r(hdl, ep, &epY);
    box bbounds(point(spX, spY), point(epX, epY));
    GEOSGeom_destroy_r(hdl, sp);
    GEOSGeom_destroy_r(hdl, ep);

    std::vector<value> results_s;
    _edge_idx->query(boost::geometry::index::contains(bbounds), std::back_inserter(results_s));
}

/**
 * Detect if a's bounding box contains b's bounding box, using the edge index.
 */
template <class T>
bool Topology::intersects(const edge* a, const T* b) {
    const GEOSGeometry* bshell = GEOSGetExteriorRing_r(hdl, b->envelope());
    GEOSGeometry* sp = GEOSGeomGetPointN_r(hdl, bshell, 0);
    GEOSGeometry* ep = GEOSGeomGetPointN_r(hdl, bshell, 3);
    double spX, spY, epX, epY;
    GEOSGeomGetX_r(hdl, sp, &spX);
    GEOSGeomGetY_r(hdl, sp, &spY);
    GEOSGeomGetX_r(hdl, ep, &epX);
    GEOSGeomGetY_r(hdl, ep, &epY);
    box bbounds(point(spX, spY), point(epX, epY));
    GEOSGeom_destroy_r(hdl, sp);
    GEOSGeom_destroy_r(hdl, ep);

    std::vector<value> results_s;
    _edge_idx->query(boost::geometry::index::intersects(bbounds), std::back_inserter(results_s));
}

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
 * Find, if it exists, the geometry from a set of geometries (others) which is the closest within
 * a specified tolerance.
 */
template <class T>
const T* closest_and_within(const GEOSGeom geom, const std::vector<T*>& others, double tolerance)
{
    const T* item = NULL;
    double previousDistance = std::numeric_limits<double>::max();
    for (T* other : others) {
        if (other && other->intersects(geom) && ST_DWithin(other->geom, geom, tolerance)) {
            double d = ST_Distance(geom, other->geom);
            if (d < previousDistance) {
                item = other;
            }
            previousDistance = d;
        }
    }
    return item;
}

/**
 * Return all items intersecting with geom in the given index.
 */
template <class IndexType, class Value>
void Topology::_intersects(IndexType* index, const GEOSGeometry* geom, std::vector<int>& edgeIds) {
    assert (index);
    assert (geom);
    assert (edgeIds.size() == 0);

    GEOSGeometry* envelope = NULL;
    const GEOSGeometry* bshell = NULL;

    std::vector<Value> results_s;

    if (GEOSGeomTypeId_r(hdl, geom) == GEOS_POINT) {
        double x, y;
        GEOSGeomGetX_r(hdl, geom, &x);
        GEOSGeomGetY_r(hdl, geom, &y);
        point pt(x, y);
        _edge_idx->query(boost::geometry::index::intersects(pt), back_inserter(results_s));
    }
    else {
        envelope = GEOSEnvelope_r(hdl, geom);
        bshell = GEOSGetExteriorRing_r(hdl, envelope);

        GEOSGeometry* sp = GEOSGeomGetPointN_r(hdl, bshell, 0);
        GEOSGeometry* ep = GEOSGeomGetPointN_r(hdl, bshell, 3);
        double spX, spY, epX, epY;
        GEOSGeomGetX_r(hdl, sp, &spX);
        GEOSGeomGetY_r(hdl, sp, &spY);
        GEOSGeomGetX_r(hdl, ep, &epX);
        GEOSGeomGetY_r(hdl, ep, &epY);
        box bbounds(point(spX, spY), point(epX, epY));
        GEOSGeom_destroy_r(hdl, sp);
        GEOSGeom_destroy_r(hdl, ep);

        _edge_idx->query(boost::geometry::index::intersects(bbounds), back_inserter(results_s));
    }

/*
    std::cout << _geos.as_string(geom) << std::endl;
    std::cout << _geos.as_string(envelope) << std::endl;
    std::cout << _geos.as_string(bshell) << std::endl;
*/

    if (envelope) {
        GEOSGeom_destroy_r(hdl, envelope);
    }

    if (results_s.size() > 0) {
        edgeIds.reserve(results_s.size()+1);
        transform(results_s.begin(), results_s.end(), back_inserter(edgeIds), [](const Value& a) {
            return a.second;
        });
    }
}

} // namespace cma

#endif // __CMA_TOPOLOGY_H
