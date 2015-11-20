#include <merge.h>

#include <memory>
#include <vector>

using namespace std;

namespace cma {

typedef vector<int> itemid_map;

void merge_topologies(Topology& t1, Topology& t2)
{
    assert (t1._transactions->empty());
    assert (t2._transactions->empty());

    unique_ptr<itemid_map> node_map(new itemid_map(t2._nodes.size(), -1));
    unique_ptr<itemid_map> edge_map(new itemid_map(t2._edges.size(), -1));
    unique_ptr<itemid_map> face_map(new itemid_map(t2._faces.size(), -1));

    // universal face stays the same even after merge
    (*face_map)[0] = 0;

    int newNodeId, nextNodeId;
    int newEdgeId, nextEdgeId;
    int newFaceId, nextFaceId;

    newNodeId = nextNodeId = t1._nodes.size();
    newEdgeId = nextEdgeId = t1._edges.size();
    newFaceId = nextFaceId = t1._faces.size();

    for (node* n : t2._nodes) {
        if (!n) continue;

        (*node_map)[n->id] = nextNodeId;
        n->id = nextNodeId++;
        t1._nodes.push_back(n);
    }

    for (edge* e : t2._edges) {
        if (!e) continue;

        (*edge_map)[e->id] = nextEdgeId;
        e->id = nextEdgeId++;
        t1._edges.push_back(e);
    }

    for (face* f : t2._faces) {
        if (!f || f->id == 0) continue;

        (*face_map)[f->id] = nextFaceId;
        f->id = nextFaceId++;
        t1._faces.push_back(f);
    }

    for (int i = newEdgeId; i < t1._edges.size(); ++i) {
        edge* e = t1._edges[i];
        if (!e) continue;

        e->start_node = (*node_map)[e->start_node];
        e->end_node   = (*node_map)[e->end_node];

        e->next_left_edge  = e->next_left_edge  < 0 ? -(*edge_map)[abs(e->next_left_edge)]  : (*edge_map)[e->next_left_edge];
        e->next_right_edge = e->next_right_edge < 0 ? -(*edge_map)[abs(e->next_right_edge)] : (*edge_map)[e->next_right_edge];

        e->abs_next_left_edge  = (*edge_map)[e->abs_next_left_edge];
        e->abs_next_right_edge = (*edge_map)[e->abs_next_right_edge];

        e->left_face  = (*face_map)[e->left_face];
        e->right_face = (*face_map)[e->right_face];
    }

    t2._empty(false);
}

}
