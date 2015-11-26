#include <merge.h>

#include <zones.h>

#include <memory>
#include <vector>
#include <boost/mpi/communicator.hpp>

using namespace std;
using namespace boost::mpi;

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

void merge_topologies(
    const vector<zone*>& zones,
    vector<Topology*>& topologies,
    vector<zone*>& new_zones,
    std::vector<Topology*>& new_topologies)
{
    communicator world;

    assert (topologies.size() % 2 == 0);
    for (int i = 0; i < topologies.size(); i+=2) {
        Topology* t1 = topologies[i];
        Topology* t2 = topologies[i+1];

        cout << "[" << world.rank() << "] merging topologies for zones #"
             << t1->zoneId() << " and #" << t2->zoneId() << endl;
        merge_topologies(*t1, *t2);

        zone* z1 = *find_if(zones.begin(), zones.end(),
            [t1](const zone* z) {
                return z->id() == t1->zoneId();
            }
        );
        zone* z2 = *find_if(zones.begin(), zones.end(),
            [t2](const zone* z) {
                return z->id() == t2->zoneId();
            }
        );

        OGREnvelope envelope = z1->envelope();
        envelope.Merge(z2->envelope());
        zone* merged_zone = new zone(t1->zoneId(), envelope);
        merged_zone->count(z1->count() + z2->count());
        new_zones.push_back(merged_zone);

        delete t2;

        cout << "[" << world.rank() << "] added new merged topology for"
             << " zone #" << t1->zoneId() << endl;

        new_topologies.push_back(t1);
    }
    topologies.clear();
}

void get_next_groups(
    vector<depth_group_t> all_groups,
    vector<depth_group_t> next_groups
)
{
    int current_depth = all_groups[0].first;
    next_groups.insert(
        begin(next_groups),
        all_groups.begin(),
        find_if_not(
            all_groups.begin(),
            all_groups.end(),
            [current_depth](const depth_group_t& g) {
                return g.first == current_depth;
            }
        )
    );

    // TODO: to speed things up, also add zones which can
    // be independently merged at other depths too
}


double width(const OGREnvelope& envelope)
{
    return envelope.MaxX - envelope.MinX;
}

double height(const OGREnvelope& envelope)
{
    return envelope.MaxY - envelope.MinY;
}

direction_type position(const OGREnvelope& e1, const OGREnvelope& e2)
{
    if (e1.MinX == e2.MinX && e1.MaxX == e2.MaxX) {
        if (e1.MaxY == e2.MinY) {
            return ABOVE;
        }
        if (e1.MinY == e2.MaxY) {
            return BELOW;
        }
    }

    if (e1.MinY == e2.MinY && e1.MaxY == e2.MaxY) {
        if (e1.MaxX == e2.MinX) {
            return RIGHT;
        }
        if (e1.MinX == e2.MaxX) {
            return LEFT;
        }
    }

    return OTHER;
}

} // namespace cma
