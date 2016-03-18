#include <zones.h>

#include <chrono>
#include <string>
#include <vector>
#include <sstream>
#include <cassert>
#include <fstream>
#include <iostream>

#include <boost/filesystem/path.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/filesystem/operations.hpp>

#include <types.h>
#include <utils.h>
#include <topology.h>

using namespace std;
using namespace boost::mpi;

namespace fs = boost::filesystem;

namespace cma {

/**
 * Number of rows and columns to start our zone finding algorithm.
 * TODO: de-hardcode
 * Note: rows*cols must equals nThreads
 */
const int rows = 2;
const int cols = 2;

void _minmax_extent(const GEOSGeometry* extent, double* minX, double* maxX,
    double* minY, double* maxY)
{
    assert (extent && GEOSGeomTypeId_r(hdl, extent) == GEOS_POLYGON);

    const GEOSCoordSequence* shell = GEOSGeom_getCoordSeq_r(
        hdl,
        GEOSGetExteriorRing_r(hdl, extent)
    );

    int numPoints =
        GEOSGeomGetNumPoints_r(hdl, GEOSGetExteriorRing_r(hdl, extent));
    assert (numPoints == 5);

    vector<double> Xs;
    vector<double> Ys;
    for (int idx = 0; idx < 4; ++idx) {
        double x, y;
        GEOSCoordSeq_getX_r(hdl, shell, idx, &x);
        GEOSCoordSeq_getY_r(hdl, shell, idx, &y);
        Xs.push_back(x);
        Ys.push_back(y);
    }

    auto minmaxX = minmax_element(Xs.begin(), Xs.end());
    auto minmaxY = minmax_element(Ys.begin(), Ys.end());

    *minX = *(minmaxX.first);
    *maxX = *(minmaxX.second);
    *minY = *(minmaxY.first);
    *maxY = *(minmaxY.second);
}

int prepare_zones(string& postgres_connect_str, GEOSHelper& geos, const GEOSGeometry* extent,
    vector<zone*>& zones, vector<depth_group_t>& grouping, int maxdepth,
    int* nextZoneId)
{
    assert (extent);
    assert (maxdepth >= 1);

    bool del_seq = false;
    if (nextZoneId == nullptr) {
        nextZoneId = new int;
        *nextZoneId = 0;
        del_seq = true;
    }

    double minX, maxX, minY, maxY;
    _minmax_extent(extent, &minX, &maxX, &minY, &maxY);

    size_t size;
    unsigned char* hex_extent =
        GEOSWKBWriter_writeHEX_r(hdl, geos.writer(), extent, &size);

    double width  = fabs(maxX - minX);
    double height = fabs(maxY - minY);

    double local_width  = fabs(width / cols);
    double local_height = fabs(height / rows);

    double sum_nb_lines = 0.;

    size_t zsize = zones.size();
    zones.resize(zsize + 4, nullptr);

    #pragma omp parallel for num_threads(4) reduction(+:sum_nb_lines)
    for (int i = 0; i < 4; ++i) {
        PG db(postgres_connect_str);
        int row, col;
        switch(i)
        {
        case 0: row = 0; col = 0; break;
        case 1: row = 0; col = 1; break;
        case 2: row = 1; col = 0; break;
        case 3: row = 1; col = 1; break;
        };
    //for (int row = 0; row < 2; ++row) {
    //    for (int col = 0; col < 2; ++col) {

        OGREnvelope local_extent;
        local_extent.MinX = minX + col * local_width;
        local_extent.MaxX = local_extent.MinX + local_width;
        local_extent.MinY = minY + row * local_height;
        local_extent.MaxY = local_extent.MinY + local_height;

        int nextId;
        nextId = (*nextZoneId)+i;
        zone* z = new zone(nextId, local_extent);

        /**
         * This will select a line count of lines within our zone.
         */
        GEOSGeometry* geom_extent = OGREnvelope2GEOSGeom(local_extent);
        string pg_geom = db.build_pg_geom(geom_extent);
        GEOSGeom_destroy_r(geos.handle(), geom_extent);

        ostringstream oss;
        oss << "SELECT COUNT(1) FROM way WHERE "
            << pg_geom << " && line2d_m AND "
            << "ST_Contains(" << pg_geom << ", line2d_m)";

        int nlines = 0;

        PGresult* res = db.query(oss.str().c_str());
        z->count(atoi(PQgetvalue(res, 0, 0)));
        PQclear(res);

        sum_nb_lines += z->count();

        zones[zsize+i] = z;
    }//}

    (*nextZoneId) += 4;

    GEOSFree_r(hdl, hex_extent);

    assert (zones.size() == 4);
    depth_group_t g = {
        maxdepth,
        {{ zones[0]->id(), zones[1]->id(), zones[2]->id(), zones[3]->id() }}
    };

    if (maxdepth == 1) {
        grouping.push_back(g);
        return sum_nb_lines;
    }

    int gIdx = 0;
    for (int zIdx = 0; zIdx < zones.size(); ++zIdx) {
        zone* z = zones[zIdx];
        if (z->count() > 15000) {
            vector<zone*> subzones;
            sum_nb_lines -= z->count();
            sum_nb_lines += prepare_zones(postgres_connect_str, geos, z->geom(), subzones, grouping, maxdepth-1, nextZoneId);

            auto oldIt = find(begin(zones), end(zones), z);
            auto newIt = zones.erase(oldIt);

            assert (g.second[gIdx] == z->id());
            g.second[gIdx] = subzones[0]->id();

            zones.insert(newIt, subzones.begin(), subzones.end());

            delete z;

            zIdx += subzones.size()-1;
        }
        ++gIdx;
    }

    grouping.push_back(g);

    if (del_seq) {
        delete nextZoneId;
    }

    return sum_nb_lines;
}

void write_zones(const std::string& filename, const vector<zone*>& zones, bool overwrite)
{
    fs::path file_path(filename);
    if (overwrite && fs::exists(file_path)) {
        fs::remove(filename);
    }

    GDALDriver* shpDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
    assert (shpDriver != NULL);
    GDALDataset *zonesDS = shpDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
    assert (zonesDS != NULL);
    OGRSpatialReference oSRS;
    oSRS.SetWellKnownGeogCS( "EPSG:3395" );
    OGRLayer* z_layer = zonesDS->CreateLayer("zones", &oSRS, wkbPolygon, NULL);
    assert (z_layer != NULL);

    OGRFieldDefn oField("zoneId", OFTInteger);
    assert (z_layer->CreateField(&oField) == OGRERR_NONE);

    for (int zIdx = 0; zIdx < zones.size(); ++zIdx) {
        OGREnvelope env = zones[zIdx]->envelope();

        OGRFeature* poFeature = OGRFeature::CreateFeature(z_layer->GetLayerDefn());
        poFeature->SetField("zoneId", zones[zIdx]->id());
        OGRPolygon* poly = new OGRPolygon();

        // per the specs, points must be added counter-clockwise
        OGRLinearRing* ring = new OGRLinearRing();
        ring->addPoint(zones[zIdx]->envelope().MinX, zones[zIdx]->envelope().MinY);
        ring->addPoint(zones[zIdx]->envelope().MaxX, zones[zIdx]->envelope().MinY);
        ring->addPoint(zones[zIdx]->envelope().MaxX, zones[zIdx]->envelope().MaxY);
        ring->addPoint(zones[zIdx]->envelope().MinX, zones[zIdx]->envelope().MaxY);
        ring->addPoint(zones[zIdx]->envelope().MinX, zones[zIdx]->envelope().MinY);    // this closes the ring

        poly->addRingDirectly(ring);

        // the feature takes ownership of the polygon, don't delete it.
        poFeature->SetGeometryDirectly(poly);

        assert (z_layer->CreateFeature(poFeature) == OGRERR_NONE);

        OGRFeature::DestroyFeature(poFeature);
    }

    GDALClose(zonesDS);
}

void register_zones(
    vector< vector<zone*> >& new_zones,
    vector<zone*>& all_zones,
    vector<zone*>& ordered_zones)
{
    for (int i = 0; i < new_zones.size(); ++i) {
        for (zone* z : new_zones[i]) {
            // insert the new zone at the original location
            // in the ordered zones
            ordered_zones.insert(
                find_if(
                    begin(ordered_zones),
                    end(ordered_zones),
                    [z](const zone* a) {
                        return z->id() == a->id();
                    }
                ),
                z
            );
        }
        all_zones.insert(
            end(all_zones),
            new_zones[i].begin(),
            new_zones[i].end()
        );
    }
    new_zones.clear();
}

GEOSGeometry* world_geom()
{
    OGREnvelope env;

    // see projected bounds from http://epsg.io/3395
    env.MinX = -20026376.39;
    env.MaxX = 20026376.39;
    env.MinY = -15496570.74;
    env.MaxY = 18764656.23;

    return OGREnvelope2GEOSGeom(env);
}

GEOSGeometry* OGREnvelope2GEOSGeom(const OGREnvelope& env)
{
    GEOSCoordSequence* corners = GEOSCoordSeq_create_r(hdl, 5, 2);
    assert (corners != NULL);

    GEOSCoordSeq_setX_r(hdl, corners, 0, env.MinX);
    GEOSCoordSeq_setY_r(hdl, corners, 0, env.MinY);

    GEOSCoordSeq_setX_r(hdl, corners, 1, env.MaxX);
    GEOSCoordSeq_setY_r(hdl, corners, 1, env.MinY);

    GEOSCoordSeq_setX_r(hdl, corners, 2, env.MaxX);
    GEOSCoordSeq_setY_r(hdl, corners, 2, env.MaxY);

    GEOSCoordSeq_setX_r(hdl, corners, 3, env.MinX);
    GEOSCoordSeq_setY_r(hdl, corners, 3, env.MaxY);

    GEOSCoordSeq_setX_r(hdl, corners, 4, env.MinX);
    GEOSCoordSeq_setY_r(hdl, corners, 4, env.MinY);

    GEOSGeometry* shell = GEOSGeom_createLinearRing_r(hdl, corners);
    assert (shell != NULL);

    GEOSGeometry* polygon = GEOSGeom_createPolygon_r(
        hdl,
        shell,
        NULL,
        0
    );

    assert (polygon != NULL);

    GEOSSetSRID_r(hdl, polygon, 3395);
    return polygon;
}

zone* get_zone_by_id(const vector<zone*>& zones, int zoneId)
{
    auto it = find_if(zones.begin(), zones.end(),
        [zoneId](const zone* z) {
            return z->id() == zoneId;
        }
    );

    if (it == zones.end()) {
        return nullptr;
    }

    return *it;
}

Topology* restore_topology(GEOSHelper* geos, zone* z, bool buildIndex)
{
    assert (z);
    assert (geos);

    ostringstream oss;
    oss << "topology-" << geos->as_hex_string(z->geom()) << ".ser";

    communicator world;

    ifstream ifs(oss.str());
    if (!ifs.is_open()) {
        cout << "[" << world.rank() << "] (file not found) error restoring topology for zone #" << z->id() << " from " << oss.str() << endl;
        return nullptr;
    }

    Topology* t = new Topology(geos);
    try {
        auto start = chrono::steady_clock::now();
        boost::archive::binary_iarchive ia(ifs);
        ia >> *t;
        auto end = chrono::steady_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "[" << world.rank() << "] restored topology for zone #" << z->id() << " from " << oss.str() << "(took: " << elapsed.count() << " ms)" << endl;
    } catch (const boost::archive::archive_exception&) {
        cout << "[" << world.rank() << "] error restoring topology for zone #" << z->id() << " from " << oss.str() << endl;
        delete t;
        t = nullptr;
    }

    ifs.close();

    if (t && buildIndex) {
        t->rebuild_indexes();
    }

    return t;
}

void save_topology(GEOSHelper* geos, zone* z, Topology* t)
{
    assert (t);
    assert (z);
    assert (geos);
    assert (z->id() == t->zoneId());

    ostringstream oss;
    oss << "topology-" << geos->as_hex_string(z->geom()) << ".ser";

    communicator world;

    cout << "saving topology for zone #" << z->id() << " to " << oss.str() << endl;

    auto start = chrono::steady_clock::now();
    ofstream ofs(oss.str());
    boost::archive::binary_oarchive oa(ofs);
    oa << *t;
    ofs.close();
    auto end = chrono::steady_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "[" << world.rank() << "] done saving zone #" << z->id() << " from " << oss.str() << "(took: " << elapsed.count() << " ms)" << endl;
}

bool restore_zones(vector<zone*>& zones, vector<depth_group_t>& groups)
{
    assert (zones.empty());
    assert (groups.empty());

    string filename = "zones.ser";

    ifstream ifs(filename);
    if (!ifs.is_open()) {
        return false;
    }

    cout << "restoring zones & groups information from " << filename << endl;

    boost::archive::binary_iarchive ia(ifs);
    ia >> zones;
    ia >> groups;

    cout << "found " << zones.size() << " zones & " << groups.size() << " groups" << endl;

    ifs.close();
    return true;
}

void save_zones(vector<zone*>& zones, vector<depth_group_t>& groups)
{
    assert (!zones.empty());
    assert (!groups.empty());

    string filename = "zones.ser";

    cout << "saving zones & groups information to " << filename << endl;

    ofstream ofs(filename);
    boost::archive::binary_oarchive oa(ofs);
    oa << zones;
    oa << groups;
    ofs.close();
}

} // namespace cma
