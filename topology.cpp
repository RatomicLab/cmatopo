#include <topology.h>

#include <cmath>
#include <cassert>
#include <algorithm>

#include <st.h>

using namespace std;

namespace cma {

Topology::Topology(GEOSContextHandle_t geos)
: _geos(geos)
, _nodes()
, _edges()
, _faces()
{
}

Topology::~Topology()
{
    delete_all(_nodes);
    delete_all(_edges);
    delete_all(_faces);
}

void Topology::TopoGeo_AddLineString(GEOSGeom line, std::vector<int>& edgeIds, double tolerance)
{
    assert (GEOSGeomTypeId(line) == GEOS_LINESTRING);

    if (tolerance < 0) {
        tolerance = _ST_MinTolerance(line);
    }

    // 1. Self-node
    GEOSGeom noded = GEOSUnaryUnion_r(_geos, line);

    // 2. Node to edges falling within tolerance distance
    vector<GEOSGeom> nearby;
    for (edge* e : _edges) {
        if (ST_DWithin(e->geom, noded, tolerance)) {
            nearby.push_back(e->geom);
        }
    }

    GEOSGeom iedges = GEOSGeom_createCollection_r(
        _geos,
        GEOS_MULTILINESTRING,
        nearby.data(),
        nearby.size()
    );

    nearby.clear();

    // 2.1 Node with existing nodes within tolerance
    for (node* n : _nodes) {
        if (ST_DWithin(n->geom, noded, tolerance)) {
            nearby.push_back(n->geom);
        }
    }

    GEOSGeom inodes = GEOSGeom_createCollection_r(
        _geos,
        GEOS_MULTIPOINT,
        nearby.data(),
        nearby.size()
    );

    nearby.clear();

    if (inodes && !GEOSisEmpty_r(_geos, inodes)) {
        noded = ST_Snap(noded, inodes, tolerance);

        int nbNodes = GEOSGetNumGeometries_r(_geos, inodes);

        for (int i = 0; i < nbNodes; ++i) {
            GEOSGeom coll = ST_Split(noded, GEOSGetGeometryN_r(_geos, inodes, i));
            assert (GEOSGetNumGeometries_r(_geos, coll) == 2);
            noded = ST_Collect(
                GEOSGeom_clone_r(_geos, GEOSGetGeometryN_r(_geos, coll, 0)),
                GEOSGeom_clone_r(_geos, GEOSGetGeometryN_r(_geos, coll, 1))
            );
        }

        noded = GEOSUnaryUnion_r(_geos, noded);
    }

    // 3. For each (now-noded) segment, insert an edge
    int nbNodes = GEOSGetNumGeometries_r(_geos, noded);
    for (int i = 0; i < nbNodes; ++i) {
        const GEOSGeometry* rec = GEOSGetGeometryN_r(_geos, noded, i);

        int start_node = TopoGeo_AddPoint(ST_StartPoint(_geos, rec), tolerance);
        int end_node = TopoGeo_AddPoint(ST_EndPoint(_geos, rec), tolerance);

        GEOSGeom sn = _nodes[start_node]->geom;
        GEOSGeom en = _nodes[end_node]->geom;

        GEOSGeom snapped = ST_SetPoint(_geos, ST_SetPoint(_geos, rec, ST_NPoints(rec)-1, en), 0, sn);

        snapped = ST_CollectionExtract(_geos, ST_MakeValid(snapped), 2);

        if (ST_IsEmpty(_geos, snapped)){
            continue;
        }

        int edgeId = -1;
        for (edge* e : _edges) {
            if (ST_Equals(e->geom, snapped)) {
                edgeId = e->id;
                break;
            }
        }

        if (edgeId == -1) {
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
    assert (geom && GEOSGeomTypeId(geom) == GEOS_LINESTRING);
    assert (GEOSisSimple_r(_geos, geom));

    edge* newEdge = new edge();
    newEdge->id = _edges.size();
    newEdge->geom = geom;
    newEdge->start_node = start_node;
    newEdge->end_node = end_node;

    bool isclosed = start_node == end_node;
    GEOSGeom cleangeom = ST_RemoveRepeatedPoints(newEdge->geom);

    GEOSGeom start_point = ST_StartPoint(_geos, cleangeom);
    GEOSGeom second_point = ST_PointN(_geos, cleangeom, 2);
    _span_t span;
    span.myaz = ST_Azimuth(_geos, start_point, second_point);
    GEOSGeom_destroy_r(_geos, start_point);
    GEOSGeom_destroy_r(_geos, second_point);

    GEOSGeom end_point = ST_EndPoint(_geos, cleangeom);
    GEOSGeom next_to_last_point = ST_PointN(_geos, cleangeom, GEOSGeomGetNumPoints_r(_geos, cleangeom)-1);
    _span_t epan;
    epan.myaz = ST_Azimuth(_geos, end_point, next_to_last_point);
    GEOSGeom_destroy_r(_geos, end_point);
    GEOSGeom_destroy_r(_geos, next_to_last_point);

    int nodes[] = {start_node, end_node};
    int i = 0;
    for (; i < 2; ++i) {
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
    // assert (ST_Equals(start_node_geom, ST_StartPoint(geom)));
    assert (end_node_geom);
    // assert (ST_Equals(end_node_geom, ST_EndPoint(geom)));

    for (node* n : _nodes) {
        assert (GEOSRelatePattern_r(_geos, n->geom, geom, "T********") == 0);
    }

    for (edge* e : _edges) {
        char* relate = GEOSRelate_r(_geos, e->geom, geom);
        assert (relate != NULL);
        if (GEOSRelatePatternMatch_r(_geos, relate, "F********") == 1) {
            GEOSFree(relate);
            continue;
        }

        assert (GEOSRelatePatternMatch_r(_geos, relate, "1FFF*FFF2") == 0);
        assert (GEOSRelatePatternMatch_r(_geos, relate, "1********") == 0);
        assert (GEOSRelatePatternMatch_r(_geos, relate, "T********") == 0);

        GEOSFree(relate);
    }

    vector<edge*> edgesToStartNode;
    vector<edge*> edgesToEndNode;

    for (edge* e : _edges) {
        edge* ne;
        GEOSGeom g = NULL;

        if (e->start_node == start_node || e->end_node == start_node) {
            ne = new edge(e);
            ne->end_node = -1;
            g = ne->geom;
            ne->geom = ST_RemoveRepeatedPoints(g);
            edgesToStartNode.push_back(ne);
        }

        if (g) {
            GEOSGeom_destroy(g);
            g = NULL;
        }

        if (e->start_node == end_node || e->end_node == end_node) {
            ne = new edge(e);
            ne->start_node = -1;
            g = ne->geom;
            ne->geom = ST_RemoveRepeatedPoints(g);
            edgesToEndNode.push_back(ne);
        }

        if (g) {
            GEOSGeom_destroy(g);
        }
    }

    if (isclosed) {
        edge* scedge = new edge(newEdge);
        scedge->start_node = -1;
        scedge->left_face = 0;
        scedge->right_face = 0;
        scedge->geom = GEOSGeom_clone_r(_geos, cleangeom);

        edge* ecedge = new edge(newEdge);
        ecedge->end_node = -1;
        ecedge->left_face = 0;
        ecedge->right_face = 0;
        ecedge->geom = GEOSGeom_clone_r(_geos, cleangeom);

        edgesToStartNode.push_back(scedge);
        edgesToEndNode.push_back(ecedge);
    }

    _find_links_to_node(start_node, edgesToStartNode, span, newEdge, isclosed);

    if (_is_null(span.nextCW)) {
        newEdge->next_right_edge = newEdge->id;
        newEdge->prev_left_edge = -newEdge->id;
    }
    else {
        newEdge->next_right_edge = span.nextCW;
        newEdge->prev_left_edge = -span.nextCCW;
    }

    _find_links_to_node(end_node, edgesToEndNode, epan, newEdge, isclosed);

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

    _edges.push_back(newEdge);

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
        // TODO: update topo
        /*
      --------------------------------------------
      -- Update topogeometries, if needed
      --------------------------------------------

        -- NOT IN THE SPECS:
        -- update TopoGeometry compositions to add newface
        sql := 'SELECT r.topogeo_id, r.layer_id FROM '
          || quote_ident(atopology)
          || '.relation r, topology.layer l '
          || ' WHERE l.topology_id = ' || topoid
          || ' AND l.level = 0 '
          || ' AND l.layer_id = r.layer_id '
          || ' AND r.element_id = ' || newedge.left_face
          || ' AND r.element_type = 3 ';
        --RAISE DEBUG 'SQL: %', sql;
        FOR rec IN EXECUTE sql
        LOOP
    #ifdef POSTGIS_TOPOLOGY_DEBUG
          RAISE DEBUG 'TopoGeometry % in layer % contained the face being split (%) - updating to contain also new face %', rec.topogeo_id, rec.layer_id, newedge.left_face, newface;
    #endif

          -- Add reference to the other face
          sql := 'INSERT INTO ' || quote_ident(atopology)
            || '.relation VALUES( ' || rec.topogeo_id
            || ',' || rec.layer_id || ',' || newface || ', 3)';
          --RAISE DEBUG 'SQL: %', sql;
          EXECUTE sql;

        END LOOP;

        */
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
            GEOSGeom_clone_r(_geos, g);
        }
        geometries.push_back(g);
    }

    GEOSGeometry* shell_geoms = GEOSGeom_createCollection_r(
        _geos,
        GEOS_MULTILINESTRING,
        geometries.data(),
        geometries.size()
    );

    GEOSGeometry* shell = GEOSGeom_createPolygon_r(
        _geos,
        shell_geoms,
        NULL,
        0
    );

    bool isccw = GEOSEqualsExact_r(_geos, shell, ST_ForceRHR(shell), 0.) == 0;

    if (faceId == 0 && !isccw) {
        return NULLint;
    }

    if (mbrOnly && faceId != 0) {
        if (isccw) {
            _faces[faceId]->mbr = GEOSEnvelope_r(_geos, shell);
        }
        return NULLint;
    }

    face* newFace = new face;
    newFace->id = _faces.size();
    if (faceId != 0 && !isccw) {
        newFace->mbr = GEOSGeom_clone_r(_geos, _faces[faceId]->mbr);
    }
    else {
        newFace->mbr = GEOSEnvelope_r(_geos, shell);
    }
    _faces.push_back(newFace);

    for (int id : newRingEdges) {
        edge* e = _edges[abs(id)];
        if (e->left_face == faceId) {
            e->left_face = newFace->id;
        }
        if (e->right_face == faceId) {
            e->right_face = newFace->id;
        }
    }

    bool ishole = (faceId != 0 && !isccw);

    vector<int> absNewRingEdges;
    transform(newRingEdges.begin(), newRingEdges.end(), absNewRingEdges.begin(), [](int id){
        return abs(id);
    });
    for (edge* e : _edges) {
        if ((e->left_face == faceId || e->right_face == faceId) &&
            !_is_in(e->id, absNewRingEdges))
        {
            GEOSGeom closestPoint = GEOSInterpolate_r(_geos, e->geom, 0.2);
            bool c = ST_Contains(shell, closestPoint);
            GEOSGeom_destroy_r(_geos, closestPoint);

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

    for (node* n : _nodes) {
        if (n->containing_face == faceId) {
            bool c = ST_Contains(shell, n->geom);
            c = ishole ? !c : c;

            if (c) {
                n->containing_face = newFace->id;
            }
        }
    }

    return newFace->id;
}


void Topology::GetRingEdges(int edgeId, vector<int>& ringEdgeIds, int maxEdges)
{
    edge* currentEdge = _edges[abs(edgeId)];

    int n = 0;
    while (true) {
        assert (currentEdge);
        ringEdgeIds.push_back(currentEdge->id);

        if (edgeId < 0) {
            if (_is_null(currentEdge->next_right_edge)) {
                break;
            }
            edgeId = currentEdge->next_right_edge;
            currentEdge = _edges[abs(edgeId)];
        }
        else {
            if (_is_null(currentEdge->next_left_edge)) {
                break;
            }
            edgeId = currentEdge->next_left_edge;
            currentEdge = _edges[abs(edgeId)];
        }

        if (!_is_null(maxEdges) && ++n > maxEdges) {
            // Max traversing limit hit.
            throw exception();
        }
    }
}

void Topology::_find_links_to_node(int nodeId, std::vector<edge*>& edges, _span_t& span, edge* newEdge, bool isclosed)
{
    int i = 0;
    for (edge* e : edges) {
        ++i;

        double az = NULLdbl;

        if (e->start_node == nodeId) {
            GEOSGeom start_point = ST_StartPoint(_geos, e->geom);
            GEOSGeom second_point = ST_PointN(_geos, e->geom, 2);
            az = ST_Azimuth(_geos, start_point, second_point);
            GEOSGeom_destroy_r(_geos, start_point);
            GEOSGeom_destroy_r(_geos, second_point);
        }
        else {
            // last edge
            GEOSGeom end_point = ST_EndPoint(_geos, e->geom);
            GEOSGeom next_to_last_point = ST_PointN(_geos, e->geom, GEOSGeomGetNumPoints_r(_geos, e->geom)-1);
            az = ST_Azimuth(_geos, end_point, next_to_last_point);
            GEOSGeom_destroy_r(_geos, end_point);
            GEOSGeom_destroy_r(_geos, next_to_last_point);

            e->id = -e->id;
        }

        assert (!_is_null(az));

        az -= span.myaz;
        if (az < 0) {
            az += 2*M_PI;
        }

        if (_is_null(span.maxaz) || az > span.maxaz) {
            span.maxaz = az;
            span.nextCCW = e->id;
            if (abs(e->id) != newEdge->id) {
                if (e->id < 0) {
                    newEdge->left_face = e->left_face;
                }
                else {
                    newEdge->left_face = e->right_face;
                }
            }
        }

        if (_is_null(span.minaz) || az < span.minaz) {
            span.minaz = az;
            span.nextCW = e->id;
            if (abs(e->id) != newEdge->id) {
                if (e->id < 0) {
                    newEdge->right_face = e->right_face;
                }
                else {
                    newEdge->right_face = e->left_face;
                }
            }
        }
    }

    span.was_isolated = isclosed ? (i < 2 ? true : false) : (i < 1 ? true : false);
}

/**
 * Mostly equivalent to the following function (topology/sql/sqlmm.sql.in):
 *   FUNCTION topology.ST_ModEdgeSplit(atopology varchar, anedge integer, apoint geometry)
 *
 * This will split an existing edge at the provided point geometry.
 */
int Topology::ST_ChangeEdgeGeom(int edgeId, const GEOSGeom acurve)
{
    assert (edgeId >= 0 && edgeId < _edges.size());
    assert (acurve != NULL);

    edge* oldEdge = _edges[edgeId];
    assert (oldEdge);

    assert (ST_Equals(ST_StartPoint(_geos, acurve), ST_StartPoint(_geos, oldEdge->geom)));

    if (oldEdge->start_node == oldEdge->end_node) {
        GEOSGeom range = ST_MakePolygon(oldEdge->geom);
        bool iscw = ST_OrderingEquals(range, ST_ForceRHR(range));

        assert (ST_NPoints(ST_RemoveRepeatedPoints(acurve)) >= 3);
        range = ST_MakePolygon(acurve);
        assert (iscw == ST_OrderingEquals(range, ST_ForceRHR(range)));
    }
    else {
        assert (ST_Equals(ST_EndPoint(_geos, acurve), ST_EndPoint(_geos, oldEdge->geom)));
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

    for (node* n : _nodes) {
        GEOSGeom coll = ST_Collect(ST_Envelope(oldEdge->geom), ST_Envelope(acurve));
        if (coll && n->id != oldEdge->start_node && n->id != oldEdge->end_node) {
            GEOSGeom ns = ST_Collect(n->geom);
            GEOSGeom r1 = NULL;
            GEOSGeom r2 = NULL;

            if (!ST_IsEmpty(_geos, ns)) {
                assert (false);
                /* implement only if necessary
                tmp1 := ST_MakeLine(ST_EndPoint(oldedge.geom), ST_StartPoint(oldedge.geom));

                rng_info.r1 := ST_MakeLine(oldedge.geom, tmp1);
                IF ST_NumPoints(rng_info.r1) < 4 THEN
                  rng_info.r1 := ST_AddPoint(rng_info.r1, ST_StartPoint(oldedge.geom));
                END IF;
                rng_info.r1 := ST_CollectionExtract(
                                   ST_MakeValid(ST_MakePolygon(rng_info.r1)), 3);

                rng_info.r2 := ST_MakeLine(acurve, tmp1);
                IF ST_NumPoints(rng_info.r2) < 4 THEN
                  rng_info.r2 := ST_AddPoint(rng_info.r2, ST_StartPoint(oldedge.geom));
                END IF;
                rng_info.r2 := ST_CollectionExtract(
                                   ST_MakeValid(ST_MakePolygon(rng_info.r2)), 3);

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
             }
        }
    }

/*

  --
  -- Check edge adjacency before
  --{

  SELECT topology._ST_AdjacentEdges(
      atopology, oldedge.start_node, anedge
    ) as pre, NULL::integer[] as post
  INTO STRICT snode_info;
#ifdef POSTGIS_TOPOLOGY_DEBUG
  RAISE DEBUG 'Bs:%', snode_info.pre;
#endif

  SELECT topology._ST_AdjacentEdges(
      atopology, oldedge.end_node, -anedge
    ) as pre, NULL::integer[] as post
  INTO STRICT enode_info;
#ifdef POSTGIS_TOPOLOGY_DEBUG
  RAISE DEBUG 'Be:%', enode_info.pre;
#endif

  --}
    */

    _edges[edgeId]->geom = acurve;

    /*
  --
  -- Check edge adjacency after
  --{

  snode_info.post := topology._ST_AdjacentEdges(
      atopology, oldedge.start_node, anedge
    );
#ifdef POSTGIS_TOPOLOGY_DEBUG
  RAISE DEBUG 'As:%', snode_info.post;
#endif

  enode_info.post := topology._ST_AdjacentEdges(
      atopology, oldedge.end_node, -anedge
    );
#ifdef POSTGIS_TOPOLOGY_DEBUG
  RAISE DEBUG 'Ae:%', enode_info.post;
#endif

  IF snode_info.pre != snode_info.post THEN
    RAISE EXCEPTION 'Edge changed disposition around start node %',
      oldedge.start_node;
  END IF;

  IF enode_info.pre != enode_info.post THEN
    RAISE EXCEPTION 'Edge changed disposition around end node %',
      oldedge.end_node;
  END IF;

  --}
    */

    if (oldEdge->left_face != 0) {
        for (face* f : _faces) {
            if (f->id == oldEdge->left_face) {
                f->mbr = ST_Envelope(ST_GetFaceGeometry(oldEdge->left_face));
            }
        }
    }
    if (oldEdge->right_face != 0 && oldEdge->right_face != oldEdge->left_face) {
        for (face* f : _faces) {
            if (f->id == oldEdge->right_face) {
                f->mbr = ST_Envelope(ST_GetFaceGeometry(oldEdge->right_face));
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
    assert (GEOSWithin(point, oldEdge->geom) == 1);

    const node* coincidentNode = get_node_at(ST_X(point), ST_Y(point));
    assert (!coincidentNode);

    node* newNode = new node();
    newNode->id = _nodes.size();
    newNode->geom = point;
    _nodes.push_back(newNode);

    GEOSGeom newedge2 = ST_Split(oldEdge->geom, point);
    GEOSGeom newedge1 = ST_GeometryN(newedge2, 0);
    newedge2 = ST_GeometryN(newedge2, 1);

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
    _edges.push_back(newEdge);

    oldEdge->geom = newedge1;
    oldEdge->next_left_edge = newEdge->id;
    oldEdge->abs_next_left_edge = newEdge->id;
    oldEdge->end_node = newNode->id;

    for (edge* e : _edges) {
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

    // TODO: update way_topo relation table

    return newNode->id;
}

void Topology::_ST_AdjacentEdges(int nodeId, int edgeId, std::vector<int>& edges)
{
    // not implemented yet
    assert(false);
    /*
    vector<int> edgestar
  ret integer[];
BEGIN
  WITH edgestar AS (
    SELECT *, count(*) over () AS cnt
    FROM topology.GetNodeEdges(atopology, anode)
  )
  SELECT ARRAY[ (
      SELECT p.edge AS prev FROM edgestar p
      WHERE p.sequence = CASE WHEN m.sequence-1 < 1 THEN cnt
                         ELSE m.sequence-1 END
    ), (
      SELECT p.edge AS prev FROM edgestar p WHERE p.sequence = ((m.sequence)%cnt)+1
    ) ]
  FROM edgestar m
  WHERE edge = anedge
  INTO ret;

  RETURN ret;
    */
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
        if (e->left_face == faceId || e->right_face == faceId) {
            geoms.push_back(e->geom);
        }
    }
    GEOSGeom coll = GEOSGeom_createCollection(GEOS_MULTILINESTRING, geoms.data(), geoms.size());
    GEOSGeom ret = ST_BuildArea(ST_Collect(coll));

    GEOSGeom_destroy(coll);

    return ret;
}


/**
 * Mostly equivalent to the following function (topology/sql/populate.sql.in):
 *   topology.TopoGeo_AddPoint
 */
int Topology::TopoGeo_AddPoint(GEOSGeom geom, double tolerance)
{
    const node* oldNode = closest_and_within(geom, _nodes, tolerance);
    if (oldNode) {
        return oldNode->id;
    }

    int id = -1;

    const edge* oldEdge = closest_and_within(geom, _edges, tolerance);
    if (oldEdge) {
        GEOSGeom point = ST_ClosestPoint(geom, oldEdge->geom);
        if (ST_Contains(oldEdge->geom, point)) {
            double snaptol = _ST_MinTolerance(point);
            GEOSGeom snapedge = ST_Snap(oldEdge->geom, point, snaptol);
            if (!ST_Equals(ST_StartPoint(_geos, oldEdge->geom), ST_StartPoint(_geos, snapedge))) {
                snapedge = ST_MakeLine(ST_StartPoint(_geos, oldEdge->geom), snapedge);
            }

            ST_ChangeEdgeGeom(oldEdge->id, snapedge);
        }

        id = ST_ModEdgeSplit(oldEdge->id, point);
    }
    else {
        // below: AddIsoNode. TODO needs review
        node* newNode = new node;
        newNode->id = _nodes.size();
        newNode->geom = geom;
        _nodes.push_back(newNode);

        id = newNode->id;

        //id := topology.ST_AddIsoNode(atopology, NULL, apoint);
    }

    assert (id != -1);
    return id;
}

/**
 * Mostly equivalent to the following function (topology/sql/sqlmm.sql.in):
 *   FUNCTION topology.ST_AddIsoNode(atopology varchar, aface integer, apoint geometry)
 */
int Topology::ST_AddIsoNode(const GEOSGeom point)
{
    assert (point != NULL);
    assert (GEOSGeomTypeId(point) == GEOS_POINT);

    // TODO
}

} // namespace cma
