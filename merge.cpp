#include <merge.h>

#include <zones.h>

#include <memory>
#include <vector>
#include <functional>
#include <boost/mpi/collectives.hpp>
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
    unique_ptr<itemid_map> relation_map(new itemid_map(t2._relations.size(), -1));

    // universal face stays the same even after merge
    (*face_map)[0] = 0;

    int nextNodeId;
    int newEdgeId, nextEdgeId;
    int nextFaceId;
    int nextTopogeoId;

    nextNodeId = t1._nodes.size();
    newEdgeId = nextEdgeId = t1._edges.size();
    nextFaceId = t1._faces.size();
    nextTopogeoId = t1._relations.size();

    for (int nodeId = 1; nodeId < t2._nodes.size(); ++nodeId) {
        node* n = t2._nodes[nodeId];
        if (n) {
            (*node_map)[n->id] = nextNodeId;
            n->id = nextNodeId;
        }
        t1._nodes.push_back(n);
        ++nextNodeId;
    }

    for (int edgeId = 1; edgeId < t2._edges.size(); ++edgeId) {
        edge* e = t2._edges[edgeId];
        if (e) {
            (*edge_map)[e->id] = nextEdgeId;
            e->id = nextEdgeId;
        }
        t1._edges.push_back(e);
        ++nextEdgeId;
    }

    for (face* f : t2._faces) {
        if (f->id == 0) continue;

        if (f) {
            (*face_map)[f->id] = nextFaceId;
            f->id = nextFaceId;
        }
        t1._faces.push_back(f);
        ++nextFaceId;
    }

    for (int topogeoId = 1; topogeoId < t2._relations.size(); ++topogeoId) {
        (*relation_map)[topogeoId] = nextTopogeoId;

        vector<relation*>* relations = t2._relations[topogeoId];
        if (relations) {
            for (relation* r : *relations) {
                r->topogeo_id = nextTopogeoId;
                switch (r->element_type)
                {
                case 2:     // LINESTRING (edge)
                    r->element_id = (*edge_map)[r->element_id];
                    break;
                case 3:     // FACE
                    r->element_id = (*face_map)[r->element_id];
                    break;
                default:
                    assert (false);
                }
            }
        }

        t1._relations.push_back(relations);
        ++nextTopogeoId;
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

int merge_topologies(
    PG& db,
    vector<zone*> zones,
    vector<Topology*>& topologies,
    vector<zone*>& new_zones,
    vector<Topology*>& new_topologies)
{
    assert (new_zones.empty());

    communicator world;

    int orphan_count = 0;

    assert (topologies.size() % 4 == 0);

    for (int i = 0; i < topologies.size()/4; ++i)
    {
        vector<zone*> temp_new_zones;
        Topology* t[2] = {nullptr, nullptr};

        for (int j = 0; j < 2; ++j)
        {
            Topology* t1 = topologies[i*4+j*2];
            Topology* t2 = topologies[i*4+j*2+1];

            orphan_count +=
                _internal_merge(db, zones, t1, t2, temp_new_zones);

            t[j] = t1;

            zones.erase(find_if(zones.begin(), zones.end(), [t1](const zone* z) {
                return z->id() == t1->zoneId();
            }));
            zones.erase(find_if(zones.begin(), zones.end(), [t2](const zone* z) {
                return z->id() == t2->zoneId();
            }));
        }
        assert (t[0] && t[1]);
        assert (temp_new_zones.size() == 2);

        // replace zones in the original vector with the new (temporary) ones
        zones.insert(zones.end(), temp_new_zones.begin(), temp_new_zones.end());

        orphan_count += _internal_merge(db, zones, t[0], t[1], new_zones);

        zones.erase(find_if(zones.begin(), zones.end(), [t](const zone* z) {
            return z->id() == t[0]->zoneId();
        }));
        zones.erase(find_if(zones.begin(), zones.end(), [t](const zone* z) {
            return z->id() == t[1]->zoneId();
        }));
        zones.push_back(new_zones[new_zones.size()-1]);

        // delete temporary zones
        for (zone* z : temp_new_zones) {
            delete z;
        }

        new_topologies.push_back(t[0]);
    }
    topologies.clear();

    int total_orphan_count;
    reduce(world, orphan_count, total_orphan_count, std::plus<int>(), 0);
    return total_orphan_count;
}

void get_next_groups(
    vector<depth_group_t>& all_groups,
    vector<depth_group_t>& next_groups
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

    all_groups.erase(
        begin(all_groups),
        find_if_not(
            begin(all_groups),
            end(all_groups),
            [current_depth](const depth_group_t& a) {
                return a.first == current_depth;
            }
        )
    );
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

int _internal_merge(
    PG& db,
    vector<zone*>& zones,
    Topology* t1,
    Topology* t2,
    vector<zone*>& new_zones)
{
    communicator world;

    merge_topologies(*t1, *t2);
    t1->rebuild_indexes();
    t1->commit();

    zone* z1 = get_zone_by_id(zones, t1->zoneId());
    zone* z2 = get_zone_by_id(zones, t2->zoneId());

    linesV orphans;
    db.get_common_lines(z1->envelope(), z2->envelope(), orphans);
    cout << "[" << world.rank() << "] adding " << orphans.size() << " lines to topology #"
         << t1->zoneId() << " (lc: " << z1->count() << "+" << z2->count() << ")" << endl;

    for (GEOSGeometry* line : orphans) {
        try {
            t1->TopoGeo_AddLineString(line, DEFAULT_TOLERANCE);
            t1->commit();
        }
        catch (const invalid_argument& ex) {
            t1->rollback();
        }
        GEOSGeom_destroy_r(hdl, line);
    }

    OGREnvelope envelope = z1->envelope();
    envelope.Merge(z2->envelope());
    zone* merged_zone = new zone(t1->zoneId(), envelope);
    merged_zone->count(z1->count() + z2->count() + orphans.size());
    new_zones.push_back(merged_zone);

    cout << "[" << world.rank() << "] added new merged topology for"
         << " zone #" << t1->zoneId() << " (lc: " << merged_zone->count() << ")" << endl;

    return orphans.size();
}

} // namespace cma
