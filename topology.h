#ifndef __CMA_TOPOLOGY_H
#define __CMA_TOPOLOGY_H

#include <limits>
#include <vector>
#include <geos_c.h>
#include <algorithm>

#include <st.h>
#include <types.h>

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
    Topology();
    ~Topology();

    template <class T>
    void delete_all(std::vector<T*>& v);

    /**
     * Add an edge (and it's endpoints) to the topology.
     */
    void TopoGeo_AddLineString(GEOSGeom line, std::vector<int>& edgeIds, double tolerance=0.);
    int ST_AddEdgeModFace(int start_node, int end_node, GEOSGeometry* geom);

    /*****************/

    int TopoGeo_AddPoint(GEOSGeom point, double tolerance);

    const node* get_node_at(double x, double y) {
        for (node* n : _nodes) {
            if (ST_X(n->geom) == x && ST_Y(n->geom) == y) {
                return n;
            }
        }
        return NULL;
    }

    int ST_ChangeEdgeGeom(int edgeId, const GEOSGeom point);
    int ST_ModEdgeSplit(int edgeId, const GEOSGeom point);
    void _ST_AdjacentEdges(int nodeId, int edgeId, std::vector<int>& edges);
    void GetNodeEdges(int nodeId, std::vector<int>& edges);
    GEOSGeom ST_GetFaceGeometry(int faceId);
    int ST_AddIsoNode(const GEOSGeom point);

private:
    std::vector<node*> _nodes;
    std::vector<edge*> _edges;
    std::vector<face*> _faces;

    template<class T>
    bool _is_in(T hay, const std::vector<T>& stack) const;

    template <class T>
    bool _is_null(T& val) const;

    int _ST_AddFaceSplit(int edgeId, int faceId, bool mbrOnly);
    void GetRingEdges(int edgeId, std::vector<int>& ringEdgeIds, int maxEdges=NULLint);
    void _find_links_to_node(int nodeId, std::vector<edge*>& edges, _span_t& span, edge* newEdge, bool isclosed);
};

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

} // namespace cma

#endif // __CMA_TOPOLOGY_H
