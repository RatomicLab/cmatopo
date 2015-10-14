#include <topology.h>

#include <cmath>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include <st.h>

using namespace std;

const int MAX_INDEX_ELEM = 100000;
const double DEFAULT_TOLERANCE = 1.;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

namespace cma {

Topology::Topology(GEOSHelper& geos)
: _nodes()
, _edges()
, _faces()
, _relations()
, _geos(geos)
, _edge_idx(new edge_idx_t)
, _edge_tol_idx(new edge_idx_t)
, _node_idx(new node_idx_t)
, _node_tol_idx(new edge_idx_t)
{
    // edges and nodes cannot have id 0
    _edges.push_back(NULL);
    _nodes.push_back(NULL);
}

Topology::~Topology()
{
    delete _edge_idx;
    delete _edge_tol_idx;
    delete _node_idx;
    delete _node_tol_idx;

    delete_all(_nodes);
    delete_all(_edges);
    delete_all(_faces);
    delete_all(_relations);
}

void Topology::TopoGeo_AddLineString(GEOSGeom line, std::vector<int>& edgeIds, double tolerance)
{
    assert (GEOSGeomTypeId_r(hdl, line) == GEOS_LINESTRING);

    if (tolerance <= 0) {
        tolerance = _ST_MinTolerance(line);
    }

    // 1. Self-node
    GEOSGeom noded = GEOSUnaryUnion_r(hdl, line);

    // 2. Node to edges falling within tolerance distance
    vector<GEOSGeom> nearby;

    vector<int> v1;
    _intersects<edge_idx_t, edge_value>(_edge_tol_idx, noded, v1);
    for (int _id : v1) {
        edge* e = _edges[_id];
        if (ST_DWithin(e->geom, noded, tolerance)) {
            nearby.push_back(e->geom);
        }
    }

    GEOSGeom iedges = GEOSGeom_createCollection_r(
        hdl,
        GEOS_MULTILINESTRING,
        nearby.data(),
        nearby.size()
    );

    nearby.clear();

    if (iedges && GEOSisEmpty_r(hdl, iedges) == 0) {
        GEOSGeometry* snapped = ST_Snap(noded, iedges, tolerance);

        noded = GEOSDifference_r(hdl, snapped, iedges);

        GEOSGeometry* set1 = GEOSIntersection_r(hdl, snapped, iedges);
        GEOSGeometry* set2 = GEOSLineMerge_r(hdl, set1);

        noded = GEOSUnion_r(hdl, noded, set2);
    }

    // 2.1 Node with existing nodes within tolerance
    v1.clear();
    _intersects<edge_idx_t, edge_value>(_node_tol_idx, noded, v1);
    for (int _id : v1) {
        node* n = _nodes[_id];
        if (ST_DWithin(n->geom, noded, tolerance)) {
            nearby.push_back(n->geom);
        }
    }

    GEOSGeom inodes = GEOSGeom_createCollection_r(
        hdl,
        GEOS_MULTIPOINT,
        nearby.data(),
        nearby.size()
    );

    nearby.clear();

    if (inodes && GEOSisEmpty_r(hdl, inodes) == 0) {
        noded = ST_Snap(noded, inodes, tolerance);

        int nbNodes = GEOSGetNumGeometries_r(hdl, inodes);

        for (int i = 0; i < nbNodes; ++i) {
            const GEOSGeometry* pt = GEOSGetGeometryN_r(hdl, inodes, i);
            assert (GEOSGeomTypeId_r(hdl, pt) == GEOS_POINT);

            /**
             * The returned collection may contain only one geometry if
             * the point (pt) is the start or end point.
             */
            noded = ST_Split(noded, pt);
        }

        noded = GEOSUnaryUnion_r(hdl, noded);
    }

    assert (noded);

    // 3. For each (now-noded) segment, insert an edge
    int nbNodes = GEOSGetNumGeometries_r(hdl, noded);
    for (int i = 0; i < nbNodes; ++i) {
        const GEOSGeometry* rec = GEOSGetGeometryN_r(hdl, noded, i);

        int start_node = TopoGeo_AddPoint(ST_StartPoint(rec), tolerance);
        int end_node = TopoGeo_AddPoint(ST_EndPoint(rec), tolerance);

        GEOSGeom sn = _nodes[start_node]->geom;
        GEOSGeom en = _nodes[end_node]->geom;

        GEOSGeom snapped = ST_SetPoint(ST_SetPoint(rec, ST_NPoints(rec)-1, en), 0, sn);
        snapped = ST_CollectionExtract(ST_MakeValid(snapped), GEOS_LINESTRING);

        if (ST_IsEmpty(snapped)){
            continue;
        }

        int edgeId = NULLint;
        for (edge* e : _edges) {
            if (!e) continue;
            if (ST_Equals(e->geom, snapped)) {
                edgeId = e->id;
                break;
            }
        }

        if (_is_null(edgeId)) {
            edgeIds.push_back(ST_AddEdgeModFace(start_node, end_node, snapped));
        }
        else {
            edgeIds.push_back(edgeId);
        }
    }
}

/**
 * Mostly equivalent to the following function (topology/sql/sqlmm.sql.in):
 *   FUNCTION topology.ST_AddEdgeModFace(atopology varchar, anode integer, anothernode integer, acurve geometry)

 * New edge created will take ownership of the GEOSGeometry.
 */
int Topology::ST_AddEdgeModFace(int start_node, int end_node, GEOSGeometry* geom)
{
    assert (!_is_null(start_node) && start_node < _nodes.size());
    assert (!_is_null(end_node) && end_node < _nodes.size());
    assert (geom && GEOSGeomTypeId_r(hdl, geom) == GEOS_LINESTRING);
    assert (GEOSisSimple_r(hdl, geom) == 1);

    edge* newEdge = new edge();
    newEdge->id = _edges.size();
    newEdge->geom = geom;
    newEdge->start_node = start_node;
    newEdge->end_node = end_node;

    bool isclosed = start_node == end_node;
    GEOSGeom cleangeom = ST_RemoveRepeatedPoints(newEdge->geom);

    GEOSGeom start_point = ST_StartPoint(cleangeom);
    GEOSGeom second_point = ST_PointN(cleangeom, 2);
    _span_t span;
    span.myaz = ST_Azimuth(start_point, second_point);
    GEOSGeom_destroy_r(hdl, start_point);
    GEOSGeom_destroy_r(hdl, second_point);

    GEOSGeom end_point = ST_EndPoint(cleangeom);
    GEOSGeom next_to_last_point = ST_PointN(cleangeom, GEOSGeomGetNumPoints_r(hdl, cleangeom)-1);
    _span_t epan;
    epan.myaz = ST_Azimuth(end_point, next_to_last_point);
    GEOSGeom_destroy_r(hdl, end_point);
    GEOSGeom_destroy_r(hdl, next_to_last_point);

    for (int i = 0; i < 2; ++i) {
        node* n = i == 0 ? _nodes[start_node] : _nodes[end_node];

        if (!_is_null(n->containing_face)) {
            if (_is_null(newEdge->left_face)) {
                newEdge->left_face = n->containing_face;
                newEdge->right_face = n->containing_face;
            }
            else {
                assert (newEdge->left_face == n->containing_face);
            }
        }
    }

    GEOSGeom start_node_geom = _nodes[start_node]->geom;
    GEOSGeom end_node_geom = _nodes[end_node]->geom;

    assert (start_node_geom);
    assert (ST_Equals(start_node_geom, ST_StartPoint(geom)));
    assert (end_node_geom);
    assert (ST_Equals(end_node_geom, ST_EndPoint(geom)));

    GEOSGeometry* geom_env = GEOSEnvelope_r(hdl, geom);

    vector<int> nodeIds;
    _intersects<node_idx_t, node_value>(_node_idx, geom_env, nodeIds);
    for (int nodeId : nodeIds) {
        node* n = _nodes[nodeId];
        char* relate = GEOSRelateBoundaryNodeRule_r(hdl, n->geom, geom, GEOSRELATE_BNR_ENDPOINT);
        assert (relate);
        assert (GEOSRelatePatternMatch_r(hdl, relate, "T********") == 0);
        GEOSFree_r(hdl, relate);
    }

    vector<int> edgeIds;
    _intersects<edge_idx_t, edge_value>(_edge_idx, geom_env, edgeIds);
    for (int edgeId : edgeIds) {
        edge* e = _edges[edgeId];

        char* relate = GEOSRelateBoundaryNodeRule_r(hdl, e->geom, geom, GEOSRELATE_BNR_ENDPOINT);
        assert (relate);
        if (GEOSRelatePatternMatch_r(hdl, relate, "F********") == 1) {
            GEOSFree_r(hdl, relate);
            continue;
        }

        assert (GEOSRelatePatternMatch_r(hdl, relate, "1FFF*FFF2") == 0);
        assert (GEOSRelatePatternMatch_r(hdl, relate, "1********") == 0);

        if (GEOSRelatePatternMatch_r(hdl, relate, "T********") == 1) {
            throw invalid_argument("SQL/MM Spatial exception - geometry crosses a node");
        }

        GEOSFree_r(hdl, relate);
    }

    GEOSGeom_destroy_r(hdl, geom_env); // double-free?

    vector<edge*> edgesToStartNode;
    vector<edge*> edgesToEndNode;

    for (edge* e : _edges) {
        if (!e) continue;

        edge* ne;
        GEOSGeom g = NULL;

        if (e->start_node == start_node || e->end_node == start_node) {
            ne = new edge(e);
            ne->start_node = e->end_node == start_node ? -1 : ne->start_node;
            ne->end_node = e->start_node == start_node ? -1 : ne->end_node;
            g = ne->geom;
            ne->geom = ST_RemoveRepeatedPoints(g);
            edgesToStartNode.push_back(ne);
        }

        if (g) {
            GEOSGeom_destroy_r(hdl, g);
            g = NULL;
        }

        if (e->start_node == end_node || e->end_node == end_node) {
            ne = new edge(e);
            ne->start_node = e->end_node == end_node ? -1 : ne->start_node;
            ne->end_node = e->start_node == end_node ? -1 : ne->end_node;
            g = ne->geom;
            ne->geom = ST_RemoveRepeatedPoints(g);
            edgesToEndNode.push_back(ne);
        }

        if (g) {
            GEOSGeom_destroy_r(hdl, g);
        }
    }

    if (isclosed) {
        edge* scedge = new edge(newEdge);
        scedge->start_node = -1;
        scedge->left_face = 0;
        scedge->right_face = 0;
        scedge->geom = GEOSGeom_clone_r(hdl, cleangeom);

        edge* ecedge = new edge(newEdge);
        ecedge->end_node = -1;
        ecedge->left_face = 0;
        ecedge->right_face = 0;
        ecedge->geom = GEOSGeom_clone_r(hdl, cleangeom);

        edgesToStartNode.push_back(scedge);
        edgesToEndNode.push_back(ecedge);
    }

    _find_links_to_node(start_node, edgesToStartNode, span, true, newEdge, isclosed);

    if (_is_null(span.nextCW)) {
        newEdge->next_right_edge = newEdge->id;
        newEdge->prev_left_edge = -newEdge->id;
    }
    else {
        newEdge->next_right_edge = span.nextCW;
        newEdge->prev_left_edge = -span.nextCCW;
    }

    _find_links_to_node(end_node, edgesToEndNode, epan, false, newEdge, isclosed);

    if (_is_null(epan.nextCW)) {
        newEdge->next_left_edge = -newEdge->id;
        newEdge->prev_right_edge = newEdge->id;
    }
    else {
        newEdge->next_left_edge = epan.nextCW;
        newEdge->prev_right_edge = -epan.nextCCW;
    }

    delete_all(edgesToStartNode);
    delete_all(edgesToEndNode);

    assert (newEdge->left_face == newEdge->right_face);
    assert (!_is_null(newEdge->left_face));

    newEdge->abs_next_left_edge = abs(newEdge->next_left_edge);
    newEdge->abs_next_right_edge = abs(newEdge->next_right_edge);
    add_edge(newEdge);

    if (abs(newEdge->prev_left_edge) != newEdge->id) {
        if (newEdge->prev_left_edge > 0) {
            edge* e = _edges[newEdge->prev_left_edge];
            e->next_left_edge = newEdge->id;
            e->abs_next_left_edge = newEdge->id;
        }
        else {
            edge* e = _edges[-newEdge->prev_left_edge];
            e->next_right_edge = newEdge->id;
            e->abs_next_right_edge = newEdge->id;
        }
    }

    if (abs(newEdge->prev_right_edge) != newEdge->id) {
        if (newEdge->prev_right_edge > 0) {
            edge* e = _edges[newEdge->prev_right_edge];
            e->next_left_edge = -newEdge->id;
            e->abs_next_left_edge = newEdge->id;
        }
        else {
            edge* e = _edges[-newEdge->prev_right_edge];
            e->next_right_edge = -newEdge->id;
            e->abs_next_right_edge = newEdge->id;
        }
    }

    if (span.was_isolated || epan.was_isolated) {
        _nodes[start_node]->containing_face = NULLint;
        _nodes[end_node]->containing_face = NULLint;
    }

    int newFaceId = _ST_AddFaceSplit(newEdge->id, newEdge->left_face, false);

    if (newFaceId == 0) {
        return newEdge->id;
    }

    if (_is_null(newFaceId)) {
        newFaceId = _ST_AddFaceSplit(-newEdge->id, newEdge->left_face, false);
    }
    else {
        _ST_AddFaceSplit(-newEdge->id, newEdge->left_face, true);
    }

    if (_is_null(newFaceId)) {
        return newEdge->id;
    }

    if (newEdge->left_face != 0) {
        for (relation* rel : _relations) {
            if (rel->element_id == newEdge->left_face && rel->element_id == 3) {
                relation* nrel = new relation();
                nrel->element_id = newFaceId;
                nrel->element_type = 3;
                _relations.push_back(nrel);
            }
        }
    }

    return newEdge->id;
}

int Topology::_ST_AddFaceSplit(int edgeId, int faceId, bool mbrOnly)
{
    if (faceId == 0 && mbrOnly) {
        return NULLint;
    }

    vector<int> newRingEdges;
    GetRingEdges(edgeId, newRingEdges);

    // IF fan.newring_edges @> ARRAY[-anedge] THEN
    if (_is_in(-edgeId, newRingEdges)) {
        return 0;
    }

    vector<GEOSGeometry*> geometries;

    for (int seq = 0; seq < newRingEdges.size(); ++seq) {
        GEOSGeom g = _edges[abs(newRingEdges[seq])]->geom;
        if (newRingEdges[seq] < 0) {
            g = ST_Reverse(g);
        }
        else {
            g = GEOSGeom_clone_r(hdl, g);
        }
        geometries.push_back(g);
    }

    GEOSGeometry* shell_geoms = GEOSGeom_createCollection_r(
        hdl,
        GEOS_MULTILINESTRING,
        geometries.data(),
        geometries.size()
    );

    GEOSGeometry* shell = GEOSPolygonize_r(hdl, geometries.data(), geometries.size());

    GEOSGeometry* forceRHR = ST_ForceRHR(shell);
    bool isccw = GEOSEqualsExact_r(hdl, shell, forceRHR, 0.) == 0;
    GEOSGeom_destroy_r(hdl, forceRHR);

    if (faceId == 0 && !isccw) {
        return NULLint;
    }

    if (mbrOnly && faceId != 0) {
        if (isccw) {
            _faces[faceId]->geom = GEOSEnvelope_r(hdl, shell);
        }

        GEOSGeom_destroy_r(hdl, shell);
        return NULLint;
    }

    face* newFace = new face;
    newFace->id = _faces.size();
    if (faceId != 0 && !isccw) {
        newFace->geom = GEOSGeom_clone_r(hdl, _faces[faceId]->geom);
    }
    else {
        newFace->geom = GEOSEnvelope_r(hdl, shell);
    }
    _faces.push_back(newFace);

    for (int id : newRingEdges) {
        edge* e = _edges[abs(id)];

        if (id > 0 && e->left_face == faceId) {
            e->left_face = newFace->id;
        }
        else if (id < 0 && e->right_face == faceId) {
            e->right_face = newFace->id;
        }
    }

    bool ishole = (faceId != 0 && !isccw);

    GEOSGeometry* shell_env = GEOSEnvelope_r(hdl, shell);

    vector<int> absNewRingEdges;
    transform(newRingEdges.begin(), newRingEdges.end(), absNewRingEdges.begin(), [](int id){
        return abs(id);
    });
    for (edge* e : _edges) {
        if (!e) continue;
        if ((e->left_face == faceId || e->right_face == faceId) &&
            !_is_in(e->id, absNewRingEdges))
        {
            GEOSGeom closestPoint = GEOSInterpolate_r(hdl, e->geom, 0.2);
            // sqlmm.sql.in:~3092
            bool c = e->intersects(shell_env) && ST_Contains(shell, closestPoint);
            GEOSGeom_destroy_r(hdl, closestPoint);

            c = ishole ? !c : c;

            if (c) {
                if (e->left_face == faceId) {
                    e->left_face = newFace->id;
                }
                if (e->right_face == faceId) {
                    e->right_face = newFace->id;
                }
            }
        }
    }

    GEOSGeom_destroy_r(hdl, shell_env);

    for (node* n : _nodes) {
        if (!n) continue;
        if (n->containing_face == faceId) {
            bool c = ST_Contains(shell, n->geom);
            c = ishole ? !c : c;

            if (c) {
                n->containing_face = newFace->id;
            }
        }
    }

    GEOSGeom_destroy_r(hdl, shell);
    return newFace->id;
}


void Topology::GetRingEdges(int edgeId, vector<int>& ringEdgeIds, int maxEdges)
{
    edge* currentEdge = _edges[abs(edgeId)];

    int n = 0;
    while (true) {
        assert (currentEdge);
        if (_is_in(edgeId, ringEdgeIds)) {
            break;
        }

        ringEdgeIds.push_back(edgeId);

        if (edgeId < 0) {
            edgeId = currentEdge->next_right_edge;
            currentEdge = _edges[abs(edgeId)];
        }
        else {
            edgeId = currentEdge->next_left_edge;
            currentEdge = _edges[abs(edgeId)];
        }

        if (!_is_null(maxEdges) && ++n > maxEdges) {
            // Max traversing limit hit.
            throw exception();
        }
    }
}

void Topology::_find_links_to_node(int nodeId, std::vector<edge*>& edges, _span_t& pan, bool span, edge* newEdge, bool isclosed)
{
    int i = 0;
    for (edge* e : edges) {
        ++i;

        double az = NULLdbl;

        if (e->start_node == nodeId) {
            GEOSGeom start_point = ST_StartPoint(e->geom);
            GEOSGeom second_point = ST_PointN(e->geom, 2);
            az = ST_Azimuth(start_point, second_point);
            GEOSGeom_destroy_r(hdl, start_point);
            GEOSGeom_destroy_r(hdl, second_point);
        }
        else {
            // last edge
            GEOSGeom end_point = ST_EndPoint(e->geom);
            GEOSGeom next_to_last_point = ST_PointN(e->geom, GEOSGeomGetNumPoints_r(hdl, e->geom)-1);
            az = ST_Azimuth(end_point, next_to_last_point);
            GEOSGeom_destroy_r(hdl, end_point);
            GEOSGeom_destroy_r(hdl, next_to_last_point);

            e->id = -e->id;
        }

        assert (!_is_null(az));

        az -= pan.myaz;
        if (az < 0.) {
            az += 2*M_PI;
        }

        if (_is_null(pan.maxaz) || az > pan.maxaz) {
            pan.maxaz = az;
            pan.nextCCW = e->id;
            if (abs(e->id) != newEdge->id) {
                if (e->id < 0) {
                    if (span) {
                        newEdge->left_face = e->left_face;
                    } else {
                        newEdge->right_face = e->left_face;
                    }
                }
                else {
                    if (span) {
                        newEdge->left_face = e->right_face;
                    } else {
                        newEdge->right_face = e->right_face;
                    }
                }
            }
        }

        if (_is_null(pan.minaz) || az < pan.minaz) {
            pan.minaz = az;
            pan.nextCW = e->id;
            if (abs(e->id) != newEdge->id) {
                if (e->id < 0) {
                    if (span) {
                        newEdge->right_face = e->right_face;
                    } else {
                        newEdge->left_face = e->right_face;
                    }
                }
                else {
                    if (span) {
                        newEdge->right_face = e->left_face;
                    } else {
                        newEdge->left_face = e->left_face;
                    }
                }
            }
        }
    }

    pan.was_isolated = isclosed ? (i < 2 ? true : false) : (i < 1 ? true : false);
}

/**
 * Mostly equivalent to the following function (topology/sql/sqlmm.sql.in):
 *   FUNCTION topology.ST_ModEdgeSplit(atopology varchar, anedge integer, apoint geometry)
 *
 * This will split an existing edge at the provided point geometry.
 */
int Topology::ST_ChangeEdgeGeom(int edgeId, const GEOSGeom acurve)
{
    assert (!_is_null(edgeId) && edgeId < _edges.size());
    assert (acurve != NULL);

    if (GEOSGeomTypeId_r(hdl, acurve) != GEOS_LINESTRING) {
        throw invalid_argument("SQL/MM Spatial exception - invalid curve");
    }

    if (GEOSisSimple_r(hdl, acurve) != 1) {
        throw invalid_argument("SQL/MM Spatial exception - curve not simple");
    }

    edge* oldEdge = _edges[edgeId];

    GEOSGeometry* sp1 = ST_StartPoint(acurve);
    GEOSGeometry* sp2 = ST_StartPoint(oldEdge->geom);
    assert (ST_Equals(sp1, sp2));
    GEOSGeom_destroy_r(hdl, sp1);
    GEOSGeom_destroy_r(hdl, sp2);

    if (oldEdge->start_node == oldEdge->end_node) {
        GEOSGeom range = ST_MakePolygon(oldEdge->geom);

        GEOSGeometry* forceRHR = ST_ForceRHR(range);
        bool iscw = ST_OrderingEquals(range, forceRHR);
        GEOSGeom_destroy_r(hdl, forceRHR);

        GEOSGeometry* g = ST_RemoveRepeatedPoints(acurve);
        assert (GEOSGeomGetNumPoints_r(hdl, g) >= 3);
        GEOSGeom_destroy_r(hdl, g);

        GEOSGeom_destroy_r(hdl, range);

        range = ST_MakePolygon(acurve);
        forceRHR = ST_ForceRHR(range);
        assert (iscw == ST_OrderingEquals(range, forceRHR));
        GEOSGeom_destroy_r(hdl, forceRHR);
    }
    else {
        GEOSGeometry* ep1 = ST_EndPoint(acurve);
        GEOSGeometry* ep2 = ST_EndPoint(oldEdge->geom);
        // cout << _geos.as_string(ep1) << " ==? " << _geos.as_string(ep2) << endl;
        if (!ST_Equals(ep1, ep2)) {
            throw invalid_argument("SQL/MM Spatial exception - end node not geometry end point.");
        }
        assert (ST_Equals(ep1, ep2));
        GEOSGeom_destroy_r(hdl, ep1);
        GEOSGeom_destroy_r(hdl, ep2);
    }

    /* below are just sanity checks, TODO later.
    string relate;
    for (node* n : _nodes) {
        if (n->id == oldEdge->start_node || n->id == oldEdge->end_node) continue;

        relate.clear();
        ST_relate(n->geom, acurve, 2, relate);
        assert (!ST_relatematch(relate, "T********"));
    }

  --
  -- h) Check if this geometry has any interaction with any existing edge
  --
  sql := 'SELECT edge_id, ST_Relate(geom,'
    || quote_literal(acurve::text)
    || '::geometry, 2) as im FROM '
    || quote_ident(atopology)
    || '.edge_data WHERE edge_id != ' || anedge || ' AND geom && '
    || quote_literal(acurve::text) || '::geometry';
  FOR rec IN EXECUTE sql LOOP -- {

    --RAISE DEBUG 'IM=%',rec.im;

    IF ST_RelateMatch(rec.im, 'F********') THEN
      CONTINUE; -- no interior-interior intersection
    END IF;

    IF ST_RelateMatch(rec.im, '1FFF*FFF2') THEN
      RAISE EXCEPTION
        'SQL/MM Spatial exception - coincident edge %', rec.edge_id;
    END IF;


    -- NOT IN THE SPECS: geometry touches an edge
    IF ST_RelateMatch(rec.im, '1********') THEN
      RAISE EXCEPTION
        'Spatial exception - geometry intersects edge %', rec.edge_id;
    END IF;

    IF ST_RelateMatch(rec.im, 'T********') THEN
      RAISE EXCEPTION
        'SQL/MM Spatial exception - geometry crosses edge %', rec.edge_id;
    END IF;

  END LOOP; -- }
    */

    vector<GEOSGeometry*> rng_info;

    GEOSGeometry* g1 = ST_Envelope(oldEdge->geom);
    GEOSGeometry* g2 = ST_Envelope(acurve);
    GEOSGeom coll = ST_Collect(g1, g2);
    GEOSGeometry* coll_env = GEOSEnvelope_r(hdl, coll);
    for (node* n : _nodes) {
        if (!n) continue;
        if (n->id != oldEdge->start_node && n->id != oldEdge->end_node && n->intersects(coll_env)) {
            rng_info.push_back(n->geom);
        }
    }
    GEOSGeom_destroy_r(hdl, coll_env);

    if (rng_info.size() > 0) {
        GEOSGeometry* nodes = GEOSGeom_createCollection_r(
            hdl,
            GEOS_MULTILINESTRING,
            rng_info.data(),
            rng_info.size()
        );

        GEOSGeometry* ep = ST_EndPoint(oldEdge->geom);
        GEOSGeometry* sp = ST_StartPoint(oldEdge->geom);
        GEOSGeometry* tmp1 = ST_MakeLine(ep, sp);
        GEOSGeom_destroy_r(hdl, ep);

        GEOSGeometry* r1 = ST_MakeLine(oldEdge->geom, tmp1);

        if (GEOSGeomGetNumPoints_r(hdl, r1) < 4) {
            GEOSGeometry* g = r1;
            r1 = ST_AddPoint(r1, sp);
            GEOSGeom_destroy_r(hdl, g);
        }

        GEOSGeometry* g = r1;
        r1 = ST_CollectionExtract(ST_MakeValid(ST_MakePolygon(r1)), GEOS_POLYGON);
        GEOSGeom_destroy_r(hdl, g);

        GEOSGeometry* r2 = ST_MakeLine(acurve, tmp1);
        if (GEOSGeomGetNumPoints_r(hdl, r2) < 4) {
            g = r2;
            r2 = ST_AddPoint(r2, sp);
            GEOSGeom_destroy_r(hdl, g);
        }

        g = r2;
        r2 = ST_CollectionExtract(ST_MakeValid(ST_MakePolygon(r2)), 3);
        GEOSGeom_destroy_r(hdl, g);

        /**
        FOR rec IN WITH
          nodes AS ( SELECT * FROM ST_Dump(rng_info.nodes) ),
          inr1 AS ( SELECT path[1] FROM nodes WHERE ST_Contains(rng_info.r1, geom) ),
          inr2 AS ( SELECT path[1] FROM nodes WHERE ST_Contains(rng_info.r2, geom) )
          ( SELECT * FROM inr1
              EXCEPT
            SELECT * FROM inr2
          ) UNION
          ( SELECT * FROM inr2
              EXCEPT
            SELECT * FROM inr1
          )
        LOOP
          RAISE EXCEPTION 'Edge motion collision at %',
                         ST_AsText(ST_GeometryN(rng_info.nodes, rec.path));
        END LOOP;
        */

        GEOSGeom_destroy_r(hdl, tmp1);
        GEOSGeom_destroy_r(hdl, sp);
    }

    GEOSGeom_destroy_r(hdl, coll);

    vector<int> preStartEdgeIds;
    vector<int> preEndEdgeIds;
    _ST_AdjacentEdges(edgeId, edgeId, preStartEdgeIds);
    _ST_AdjacentEdges(edgeId, edgeId, preEndEdgeIds);

    GEOSGeometry* g = _edges[edgeId]->geom;
    _edges[edgeId]->geom = acurve;
    if (g) {
        GEOSGeom_destroy_r(hdl, g);
    }

    vector<int> postStartEdgeIds;
    vector<int> postEndEdgeIds;
    _ST_AdjacentEdges(edgeId, edgeId, postStartEdgeIds);
    _ST_AdjacentEdges(edgeId, edgeId, postEndEdgeIds);

    assert (preStartEdgeIds == postStartEdgeIds);
    assert (preEndEdgeIds == postEndEdgeIds);

    if (oldEdge->left_face != 0) {
        for (face* f : _faces) {
            if (f->id == oldEdge->left_face) {
                g = f->geom;
                f->geom = ST_Envelope(ST_GetFaceGeometry(oldEdge->left_face));
                if (g) {
                    GEOSGeom_destroy_r(hdl, g);
                }
            }
        }
    }
    if (oldEdge->right_face != 0 && oldEdge->right_face != oldEdge->left_face) {
        for (face* f : _faces) {
            if (f->id == oldEdge->right_face) {
                g = f->geom;
                f->geom = ST_Envelope(ST_GetFaceGeometry(oldEdge->right_face));
                if (g) {
                    GEOSGeom_destroy_r(hdl, g);
                }
            }
        }
    }

    return edgeId;
}

/**
 * Mostly equivalent to the following function (topology/sql/sqlmm.sql.in):
 *   FUNCTION topology.ST_ModEdgeSplit(atopology varchar, anedge integer, apoint geometry)
 */
int Topology::ST_ModEdgeSplit(int edgeId, const GEOSGeom point)
{
    assert (edgeId >= 0 && edgeId < _edges.size());
    assert (point);

    edge* oldEdge = _edges[edgeId];
    assert (oldEdge);

    if (GEOSWithin_r(hdl, point, oldEdge->geom) != 1) {
        throw invalid_argument("SQL/MM Spatial exception - point not on edge");
    }

    const node* coincidentNode = get_node_at(ST_X(point), ST_Y(point));
    assert (!coincidentNode);

    node* newNode = new node();
    newNode->id = _nodes.size();
    newNode->geom = point;
    add_node(newNode);

    GEOSGeometry* tmp;

    GEOSGeom newedge2 = ST_Split(oldEdge->geom, point);
    GEOSGeom newedge1 = ST_GeometryN(newedge2, 0);

    tmp = newedge2;
    newedge2 = ST_GeometryN(newedge2, 1);
    GEOSGeom_destroy_r(hdl, tmp);

    assert (newedge1 && newedge2);

    edge* newEdge = new edge();
    newEdge->id = _edges.size();
    newEdge->start_node = newNode->id;
    newEdge->end_node = oldEdge->end_node;
    newEdge->next_left_edge = oldEdge->next_left_edge == -edgeId ? -newEdge->id : oldEdge->next_left_edge;
    newEdge->next_right_edge = -edgeId;
    newEdge->left_face = oldEdge->left_face;
    newEdge->right_face = oldEdge->right_face;
    newEdge->geom = newedge2;

    newEdge->abs_next_left_edge = abs(newEdge->next_left_edge);
    newEdge->abs_next_right_edge = abs(newEdge->next_right_edge);
    add_edge(newEdge);

    tmp = oldEdge->geom;
    oldEdge->geom = newedge1;
    GEOSGeom_destroy_r(hdl, tmp);
    oldEdge->next_left_edge = newEdge->id;
    oldEdge->abs_next_left_edge = newEdge->id;
    oldEdge->end_node = newNode->id;

    for (edge* e : _edges) {
        if (!e) continue;
        if (e->id != newEdge->id) {
            if (e->next_right_edge == -edgeId && e->start_node == oldEdge->end_node) {
                e->next_right_edge = -newEdge->id;
                e->abs_next_right_edge = newEdge->id;
            }
            if (e->next_left_edge == -edgeId && e->end_node == oldEdge->end_node) {
                e->next_left_edge = -newEdge->id;
                e->abs_next_left_edge = newEdge->id;
            }
        }
    }

    for (relation* r : _relations) {
        if (abs(r->element_id) == edgeId && r->element_type == 2) {
            relation* nrel = new relation;
            nrel->element_id = r->element_id < 0 ? -newEdge->id : newEdge->id;
            nrel->element_type = r->element_type;
            _relations.push_back(nrel);
        }
    }

    return newNode->id;
}

void Topology::_ST_AdjacentEdges(int nodeId, int edgeId, vector<int>& edgeIds)
{
    vector<int> edgeStar;
    GetNodeEdges(nodeId, edgeStar);

    int msequence = 0;
    int psequence = 0;

    for (int mid : edgeStar) {
        ++msequence;

        if (mid == edgeId) {
            for (int pid : edgeStar) {
                ++psequence;

                if (psequence == (msequence-1<1 ? edgeStar.size() : msequence-1) ||
                    psequence == (msequence%edgeStar.size())+1)
                {
                    edgeIds.push_back(pid);
                }
            }
        }
    }
}

void Topology::GetNodeEdges(int nodeId, vector<int>& edgeIds)
{
    typedef pair<double, int> az_edge_pair;
    vector< az_edge_pair > azs;

    for (edge* e : _edges) {
        if (!e) continue;
        if (e->start_node == nodeId || e->end_node == nodeId) {
            GEOSGeometry* geom = ST_RemoveRepeatedPoints(e->geom);

            double az;
            if (e->start_node == nodeId) {
                GEOSGeometry* p1 = ST_StartPoint(geom);
                GEOSGeometry* p2 = ST_PointN(geom, 2);
                az = ST_Azimuth(p1, p2);
                GEOSGeom_destroy_r(hdl, p1);
                GEOSGeom_destroy_r(hdl, p2);

                azs.push_back(make_pair(az, e->id));
            }
            else {
                GEOSGeometry* p1 = ST_EndPoint(geom);
                GEOSGeometry* p2 = ST_PointN(geom, GEOSGeomGetNumPoints_r(hdl, geom)-1);
                az = ST_Azimuth(p1, p2);
                GEOSGeom_destroy_r(hdl, p1);
                GEOSGeom_destroy_r(hdl, p2);

                azs.push_back(make_pair(az, -e->id));
            }

            GEOSGeom_destroy_r(hdl, geom);
        }
    }

    // sort the edges by azimuth values
    sort(azs.begin(), azs.end(), [](const az_edge_pair& a, const az_edge_pair& b) {
        return a.first < b.first;
    });

    transform(azs.begin(), azs.end(), std::back_inserter(edgeIds), [](const az_edge_pair& a) {
        return a.second;
    });
}

/**
 * Mostly equivalent to the following function (topology/sql/sqlmm.sql.in):
 *   FUNCTION topology.ST_ModEdgeSplit(atopology varchar, anedge integer, apoint geometry)
 */
GEOSGeom Topology::ST_GetFaceGeometry(int faceId)
{
    // face 0 is invalid (universal face)
    assert (faceId > 0 && faceId < _faces.size());

    vector<GEOSGeom> geoms;
    for (edge* e : _edges) {
        if (!e) continue;
        if (e->left_face == faceId || e->right_face == faceId) {
            geoms.push_back(e->geom);
        }
    }

    GEOSGeom coll = GEOSGeom_createCollection_r(
        hdl,
        GEOS_MULTILINESTRING,
        geoms.data(),
        geoms.size()
    );

    GEOSGeom ret = ST_BuildArea(coll);
    GEOSGeom_destroy_r(hdl, coll);
    return ret;
}

/**
 * Mostly equivalent to the following function (topology/sql/populate.sql.in):
 *   topology.TopoGeo_AddPoint
 */
int Topology::TopoGeo_AddPoint(GEOSGeom geom, double tolerance)
{
    assert (geom);
    assert (GEOSGeomTypeId_r(hdl, geom) == GEOS_POINT);

    if (tolerance == 0.) {
        tolerance = _ST_MinTolerance(geom);
    }

    const node* oldNode = closest_and_within_node(geom, tolerance);
    if (oldNode) {
        return oldNode->id;
    }

    int id = -1;

    const edge* oldEdge = closest_and_within_edge(geom, tolerance);
    if (oldEdge) {
        GEOSGeom point = ST_ClosestPoint(oldEdge->geom, geom);

        if (!ST_Contains(oldEdge->geom, point)) {
            double snaptol = _ST_MinTolerance(point);

            GEOSGeom snapedge = ST_Snap(oldEdge->geom, point, snaptol);
            GEOSGeom sp = ST_StartPoint(oldEdge->geom);
            GEOSGeom ep = ST_StartPoint(snapedge);

            if (!ST_Equals(sp, ep)) {
                GEOSGeom_destroy_r(hdl, sp);
                sp = ST_StartPoint(oldEdge->geom);
                snapedge = ST_MakeLine(sp, snapedge);
            }

            // snapedge must not be destroyed since it gets assigned to oldEdge->geom
            ST_ChangeEdgeGeom(oldEdge->id, snapedge);

            GEOSGeom_destroy_r(hdl, sp);
            GEOSGeom_destroy_r(hdl, ep);
        }

        // point must not be destroyed as it gets assigned to a new node's geometry
        id = ST_ModEdgeSplit(oldEdge->id, point);
    }
    else {
        id = ST_AddIsoNode(NULLint, geom);
    }

    assert (id != -1);
    return id;
}

/**
 * Mostly equivalent to the following function (topology/sql/sqlmm.sql.in):
 *   FUNCTION topology.ST_AddIsoNode(atopology varchar, aface integer, apoint geometry)
 * Most code pertaining to faces was not ported as we don't need it (yet).
 */
int Topology::ST_AddIsoNode(int faceId, const GEOSGeom point)
{
    assert (point != NULL);
    assert (GEOSGeomTypeId_r(hdl, point) == GEOS_POINT);

    for (node* n : _nodes) {
        if (!n) continue;
        if (ST_Equals(n->geom, point)) {
            throw invalid_argument("SQL/MM Spatial exception - coincident node");
        }
    }

    vector<int> edgeIds;
    _intersects<edge_idx_t, edge_value>(_edge_idx, point, edgeIds);
    if (any_of(edgeIds.begin(), edgeIds.end(), [this, point](int edgeId) {
        return GEOSIntersects_r(hdl, this->_edges[edgeId]->geom, point) == 1;
    }))
    {
        throw invalid_argument("SQL/MM Spatial exception - edge crosses node.");
    }

    int containing_face = NULLint;
    for (face* f : _faces) {
        if (!f) continue;
        if (!_is_null(faceId) && faceId != 0 && f->id != faceId) continue;

        if (f->id > 0 && f->intersects(point) && ST_Contains(ST_GetFaceGeometry(f->id), point)) {
            containing_face = f->id;
            break;
        }
    }

    if (!_is_null(faceId)) {
        if (faceId == 0) {
            if (!_is_null(containing_face)) {
                throw invalid_argument("SQL/MM Spatial exception - within face (not universe)");
            }
            else {
                containing_face = 0;
            }
        }
        else {
            if (_is_null(containing_face) || containing_face != faceId) {
                throw invalid_argument("SQL/MM Spatial exception - not within face");
            }
        }
    }
    else {
        containing_face = _is_null(containing_face) ? 0 : containing_face;
    }

    node* newNode = new node;
    newNode->id = _nodes.size();
    newNode->geom = point;
    newNode->containing_face = containing_face;
    add_node(newNode);

    return newNode->id;
}

void Topology::output_nodes() const
{
    cout << "node count: " << _nodes.size()-1 << endl;
    cout << "node_id | containing_face | geom (as text)" << endl;
    for (const node* n : _nodes) {
        if (!n) continue;
        cout << n->id << " | " << n->containing_face << " | " << _geos.as_string(n->geom) << endl;
    }
}

void Topology::add_edge(edge* e)
{
    assert (e);
    assert (e->id == _edges.size());
    assert (e->geom);
    assert (GEOSGeomTypeId_r(hdl, e->geom) == GEOS_LINESTRING);

    bg::model::linestring<point> edgeLS;

    double x, y;
    for (int i = 0; i < GEOSGeomGetNumPoints_r(hdl, e->geom); ++i) {
        GEOSGeometry* pt = GEOSGeomGetPointN_r(hdl, e->geom, i);

        GEOSGeomGetX_r(hdl, pt, &x);
        GEOSGeomGetY_r(hdl, pt, &y);
        edgeLS.push_back(point(x, y));

        GEOSGeom_destroy_r(hdl, pt);
    }

    bg::model::box<point> bounds;
    bg::envelope(edgeLS, bounds);

    bg::model::box<point> tolbounds;
    bg::envelope(edgeLS, tolbounds);
    point& pt1 = tolbounds.min_corner();
    point& pt2 = tolbounds.max_corner();
    bg::set<0>(pt1, bg::get<0>(pt1)-DEFAULT_TOLERANCE/2);
    bg::set<1>(pt1, bg::get<1>(pt1)-DEFAULT_TOLERANCE/2);
    bg::set<0>(pt2, bg::get<0>(pt2)+DEFAULT_TOLERANCE/2);
    bg::set<1>(pt2, bg::get<1>(pt2)+DEFAULT_TOLERANCE/2);

    _edges.push_back(e);
    _edge_idx->insert(make_pair(bounds, e->id));
    _edge_tol_idx->insert(make_pair(tolbounds, e->id));
}

void Topology::add_node(node* n)
{
    assert (n);
    assert (n->geom);
    assert (n->id == _nodes.size());
    assert (GEOSGeomTypeId_r(hdl, n->geom) == GEOS_POINT);

    double x, y;
    GEOSGeomGetX_r(hdl, n->geom, &x);
    GEOSGeomGetY_r(hdl, n->geom, &y);
    point pt(x, y);

    // Index the buffer's envelope around the node.
    // See: www.boost.org/doc/libs/1_59_0/libs/geometry/doc/html/geometry/reference/strategies/strategy_buffer_point_square.html
    bg::strategy::buffer::point_square point_strategy;
    bg::strategy::buffer::distance_symmetric<double> distance_strategy(DEFAULT_TOLERANCE/2);    // radius is half the tolerance
    bg::strategy::buffer::join_round join_strategy;
    bg::strategy::buffer::end_round end_strategy;
    bg::strategy::buffer::side_straight side_strategy;

    bg::model::multi_polygon<polygon> result;
    bg::buffer(pt, result,
               distance_strategy, side_strategy,
               join_strategy, end_strategy, point_strategy);
    box tolbounds;
    bg::envelope(result, tolbounds);

    _nodes.push_back(n);
    _node_idx->insert(make_pair(pt, n->id));
    _node_tol_idx->insert(make_pair(tolbounds, n->id));
}

void Topology::output_edges() const
{
    cout << "edge count: " << _edges.size()-1 << endl;
    cout << "edge_id | start_node | end_node | next_left_edge | abs_next_left_edge | next_right_edge | abs_next_right_edge | left_face | right_face | geom (as text)" << endl;
    for (const edge* e : _edges) {
        if (!e) continue;
        cout << e->id << " | "
             << e->start_node << " | "
             << e->end_node << " | "
             << e->next_left_edge << " | "
             << e->abs_next_left_edge << " | "
             << e->next_right_edge << " | "
             << e->abs_next_right_edge << " | "
             << e->left_face << " | "
             << e->right_face << " | "
             << _geos.as_string(e->geom)
             << endl;
    }
}

/**
 * Find, if it exists, the node which is the closest to geom within a specified tolerance.
 */
const node* Topology::closest_and_within_node(const GEOSGeometry* geom, double tolerance)
{
    assert (geom);
    assert (GEOSGeomTypeId_r(hdl, geom) == GEOS_POINT);

    double x, y;
    GEOSGeomGetX_r(hdl, geom, &x);
    GEOSGeomGetY_r(hdl, geom, &y);

    vector<node_value> results_s;
    _node_idx->query(bgi::nearest(point(x, y), 1), back_inserter(results_s));

    vector<int> nodeIds;
    if (results_s.size() > 0) {
        nodeIds.reserve(results_s.size()+1);
        transform(results_s.begin(), results_s.end(), back_inserter(nodeIds), [](const node_value& a) {
            return a.second;
        });
    }

    // remove nodes not within tolerance
    nodeIds.erase(remove_if(nodeIds.begin(), nodeIds.end(), [this, geom, tolerance](int nodeId) {
        return !ST_DWithin(this->_nodes[nodeId]->geom, geom, tolerance);
    }), nodeIds.end());

    if (nodeIds.size() == 0) {
        return NULL;
    }

    if (nodeIds.size() == 1) {
        return _nodes[nodeIds[0]];
    }

    // compute distance for remaining nodes
    vector<double> distances(nodeIds.size()+1);
    transform(nodeIds.begin(), nodeIds.end(), back_inserter(distances), [this, geom](int nodeId) {
        return ST_Distance(geom, this->_nodes[nodeId]->geom);
    });

    // return the node with the smallest distance to geom
    return _nodes[nodeIds[distance(distances.begin(), min_element(distances.begin(), distances.end()))]];
}

/**
 * Find, if it exists, the edge which is the closest to geom within a specified tolerance.
 */
const edge* Topology::closest_and_within_edge(const GEOSGeometry* geom, double tolerance)
{
    vector<int> edgeIds;
    _intersects<edge_idx_t, edge_value>(_edge_idx, geom, edgeIds);

    // remove edges not within tolerance
    edgeIds.erase(remove_if(edgeIds.begin(), edgeIds.end(), [this, geom, tolerance](int edgeId) {
        return !ST_DWithin(this->_edges[edgeId]->geom, geom, tolerance);
    }), edgeIds.end());

    if (edgeIds.size() == 0) {
        return NULL;
    }

    if (edgeIds.size() == 1) {
        return _edges[edgeIds[0]];
    }

    // compute distance for remaining edges
    vector<double> distances(edgeIds.size()+1);
    transform(edgeIds.begin(), edgeIds.end(), back_inserter(distances), [this, geom](int edgeId) {
        return ST_Distance(geom, this->_edges[edgeId]->geom);
    });

    // return the edge with the smallest distance to geom
    return _edges[edgeIds[distance(distances.begin(), min_element(distances.begin(), distances.end()))]];
}

void GEOM2BOOSTMLS(const GEOSGeometry* in, multi_linestring& mls)
{
    assert (in);
    assert (GEOSGeomTypeId_r(hdl, in) == GEOS_MULTILINESTRING);
    assert (mls.size() == 0);

    double x, y;
    mls.resize(GEOSGetNumGeometries_r(hdl, in));

    for (int l = 0; l < GEOSGetNumGeometries_r(hdl, in); ++l)
    {
        const GEOSGeometry* line = GEOSGetGeometryN_r(hdl, in, l);

        for (int p = 0; p < GEOSGeomGetNumPoints_r(hdl, line); ++p)
        {
            GEOSGeometry* pt = GEOSGeomGetPointN_r(hdl, line, p);
            GEOSGeomGetX_r(hdl, pt, &x);
            GEOSGeomGetY_r(hdl, pt, &y);
            GEOSGeom_destroy_r(hdl, pt);

            boost::geometry::append(mls[l], point(x, y));
        }
    }
}

void GEOM2BOOSTLS(const GEOSGeometry* in, linestring& ls)
 {
     assert (in);
     assert (GEOSGeomTypeId_r(hdl, in) == GEOS_LINESTRING || GEOSGeomTypeId_r(hdl, in) == GEOS_LINEARRING);
     assert (ls.size() == 0);

     double x, y;
     for (int i = 0; i < GEOSGeomGetNumPoints_r(hdl, in); ++i) {
         GEOSGeometry* pt = GEOSGeomGetPointN_r(hdl, in, i);

         GEOSGeomGetX_r(hdl, pt, &x);
         GEOSGeomGetY_r(hdl, pt, &y);
         ls.push_back(point(x, y));

         GEOSGeom_destroy_r(hdl, pt);
     }
 }

} // namespace cma
