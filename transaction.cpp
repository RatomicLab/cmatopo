#include <transaction.h>

#include <topology.h>

using namespace std;

namespace cma {

EdgeTransaction::EdgeTransaction(Topology& topology, edge* e)
: TopologyTransaction(topology)
, _edge(new edge(e, false))
{
}

EdgeTransaction::~EdgeTransaction()
{
    if (_topology._edges[_edge->id]->geom == _edge->geom
     || _is_in(_edge->geom, *_topology._tr_track_geom))
    {
        _edge->geom = nullptr;
    }
    else {
        _topology._tr_track_geom->insert(_edge->geom);
    }
    delete _edge;
}

void EdgeTransaction::rollback()
{
    edge* oldEdge = _topology._edges[_edge->id];
    _topology._edges[_edge->id] = _edge;
    _edge = oldEdge;
}

NodeTransaction::NodeTransaction(Topology& topology, node* n)
: TopologyTransaction(topology)
, _node(new node(n, false))
{
}

NodeTransaction::~NodeTransaction()
{
    if (_topology._nodes[_node->id]->geom == _node->geom
     || _is_in(_node->geom, *_topology._tr_track_geom))
    {
        _node->geom = nullptr;
    }
    else {
        _topology._tr_track_geom->insert(_node->geom);
    }
    delete _node;
}

void NodeTransaction::rollback()
{
    node* oldNode = _topology._nodes[_node->id];
    _topology._nodes[_node->id] = _node;
    _node = oldNode;
}

FaceTransaction::FaceTransaction(Topology& topology, face* f)
: TopologyTransaction(topology)
, _face(new face(f, false))
{
}

FaceTransaction::~FaceTransaction()
{
    if (_topology._faces[_face->id]->geom == _face->geom
     || _is_in(_face->geom, *_topology._tr_track_geom))
    {
        _face->geom = nullptr;
    }
    else {
        _topology._tr_track_geom->insert(_face->geom);
    }
    delete _face;
}

void FaceTransaction::rollback()
{
    face* oldFace = _topology._faces[_face->id];
    _topology._faces[_face->id] = _face;
    _face = oldFace;
}

AddFaceIndexTransaction::AddFaceIndexTransaction(
    Topology& topology,
    vector<edgeid_set_ptr>* index,
    int faceId,
    int edgeId
)
: TopologyTransaction(topology)
, _index(index)
, _faceId(faceId)
, _edgeId(edgeId)
{
}

AddFaceIndexTransaction::~AddFaceIndexTransaction()
{
}

void AddFaceIndexTransaction::rollback()
{
    size_t nelem = (*_index)[_faceId]->erase(_edgeId);
    assert (nelem == 1);

    // invalidate face geometry cache
    if ((*_topology._face_geometries)[_faceId]) {
        GEOSGeom_destroy_r(hdl, (*_topology._face_geometries)[_faceId]);
        (*_topology._face_geometries)[_faceId] = nullptr;
    }
}

RemoveFaceIndexTransaction::RemoveFaceIndexTransaction(
    Topology& topology,
    vector<edgeid_set_ptr>* index,
    int faceId,
    int edgeId
)
: TopologyTransaction(topology)
, _index(index)
, _faceId(faceId)
, _edgeId(edgeId)
{
}

RemoveFaceIndexTransaction::~RemoveFaceIndexTransaction()
{
}

void RemoveFaceIndexTransaction::rollback()
{
    (*_index)[_faceId]->insert(_edgeId);

    // invalidate face geometry cache
    if ((*_topology._face_geometries)[_faceId]) {
        GEOSGeom_destroy_r(hdl, (*_topology._face_geometries)[_faceId]);
        (*_topology._face_geometries)[_faceId] = nullptr;
    }
}

AddRelationTransaction::AddRelationTransaction(
    Topology& topology,
    int topogeoId,
    int elementId,
    int elementType)
: TopologyTransaction(topology)
, _topogeoId(topogeoId)
, _elementId(elementId)
, _elementType(elementType)
{
}

AddRelationTransaction::~AddRelationTransaction()
{
}

void AddRelationTransaction::rollback()
{
    assert (_topogeoId < _topology._relations.size());

    vector<relation*>* relations = _topology._relations[_topogeoId];

    auto it = find_if(
        relations->begin(),
        relations->end(),
        [this](const relation* r) {
            return r->element_id == this->_elementId &&
                r->element_type == this->_elementType;
        }
    );
    assert (it != relations->end());

    delete *it;
    relations->erase(it);

    if (relations->empty()) {
        delete relations;
        _topology._relations[_topogeoId] = nullptr;
    }
}

} // namespace cma
