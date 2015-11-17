#include <transaction.h>

#include <topology.h>

using namespace std;

namespace cma {

EdgeTransaction::EdgeTransaction(Topology& topology, edge* e)
: TopologyTransaction(topology)
, _edge(new edge(e))
{
}

EdgeTransaction::~EdgeTransaction()
{
    delete _edge;
}

void EdgeTransaction::commit()
{
    // don't destroy _edge->geom, it will be destroyed when
    // _edge is deleted.
}

void EdgeTransaction::rollback()
{
    edge* oldEdge = _topology._edges[_edge->id];
    _topology._edges[_edge->id] = _edge;
    _edge = oldEdge;
}

NodeTransaction::NodeTransaction(Topology& topology, node* n)
: TopologyTransaction(topology)
, _node(new node(n))
{
}

NodeTransaction::~NodeTransaction()
{
    delete _node;
}

void NodeTransaction::commit()
{
    // don't destroy _node->geom, it will be destroyed when
    // _node is deleted.
}

void NodeTransaction::rollback()
{
    node* oldNode = _topology._nodes[_node->id];
    _topology._nodes[_node->id] = _node;
    _node = oldNode;
}

FaceTransaction::FaceTransaction(Topology& topology, face* f)
: TopologyTransaction(topology)
, _face(new face(f))
{
}

FaceTransaction::~FaceTransaction()
{
    delete _face;
}

void FaceTransaction::commit()
{
    // don't destroy _face->geom, it will be destroyed when
    // _face is deleted.
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
    // cout << "AddFaceIndexTransaction::rollback(): removing edge " << _edgeId << " from face " << _faceId << " index" << endl;
    size_t nelem = (*_index)[_faceId]->erase(_edgeId);
    assert (nelem == 1);
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
    // cout << "RemoveFaceIndexTransaction::rollback(): re-adding edge " << _edgeId << " to face " << _faceId << " index" << endl;
    (*_index)[_faceId]->insert(_edgeId);
}

AddRelationTransaction::AddRelationTransaction(
    Topology& topology,
    int relationId
)
: TopologyTransaction(topology)
, _relationId(relationId)
{
}

AddRelationTransaction::~AddRelationTransaction()
{
}

void AddRelationTransaction::rollback()
{
    assert (_relationId < _topology._relations.size());

    relation* rel = _topology._relations[_relationId];
    assert (rel && rel->id == _relationId);
    _topology._relations[_relationId] = nullptr;
    delete rel;
}

} // namespace cma
