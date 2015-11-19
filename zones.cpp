#include <zones.h>

#include <omp.h>
#include <assert.h>

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <utils.h>

using namespace std;

namespace fs = boost::filesystem;

namespace cma {

/**
 * Number of rows and columns to start our zone finding algorithm.
 * TODO: de-hardcode
 * Note: rows*cols must equals nThreads
 */
const int rows = 2;
const int cols = 2;

void _minmax_extent(const GEOSGeometry* extent, double* minX, double* maxX, double* minY, double* maxY)
{
    assert (extent && GEOSGeomTypeId_r(hdl, extent) == GEOS_POLYGON);

    const GEOSCoordSequence* shell = GEOSGeom_getCoordSeq_r(
        hdl,
        GEOSGetExteriorRing_r(hdl, extent)
    );

    int numPoints = GEOSGeomGetNumPoints_r(hdl, GEOSGetExteriorRing_r(hdl, extent));
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

// void prepare_zones(const linesV& lines, const OGREnvelope& extent, vector<zoneInfo*>& zones, int maxdepth)
void prepare_zones(GEOSHelper& geos, const GEOSGeometry* extent, vector<zoneInfo*>& zones, int maxdepth)
{
    assert (extent);
    assert (maxdepth >= 1);
    assert (!omp_in_parallel());

    double minX, maxX, minY, maxY;
    _minmax_extent(extent, &minX, &maxX, &minY, &maxY);

    size_t size;
    unsigned char* hex_extent = GEOSWKBWriter_writeHEX_r(hdl, geos.writer(), extent, &size);

    omp_set_num_threads(rows*cols);
    int nThreads = get_nb_threads();
    double width  = fabs(maxX - minX);
    double height = fabs(maxY - minY);

#if 0
    cout <<
        "MinX: " << minX << " "
        "MaxX: " << maxX << " "
        "MinY: " << minY << " "
        "MaxY: " << maxY <<
        endl;
#endif

    assert (rows*cols == nThreads);
    double local_width  = fabs(width / cols);
    double local_height = fabs(height / rows);

    double sum_nb_lines = 0.;

    #pragma omp parallel reduction(+:sum_nb_lines)
    {
        PG db("postgresql://postgres@localhost/postgis");

        int row = int(floor(omp_get_thread_num() / cols));
        int col = omp_get_thread_num() % cols;

        OGREnvelope local_extent;
        local_extent.MinX = minX + col * local_width;
        local_extent.MaxX = local_extent.MinX + local_width;
        local_extent.MinY = minY + row * local_height;
        local_extent.MaxY = local_extent.MinY + local_height;

        /*
        ostringstream local_wkt_extent;
        local_wkt_extent << "LINESTRING (";
        local_wkt_extent << local_extent.MinX << " " << local_extent.MinY << ", ";  // lower left
        local_wkt_extent << local_extent.MaxX << " " << local_extent.MinY << ", ";  // lower right
        local_wkt_extent << local_extent.MaxX << " " << local_extent.MaxY << ", ";  // upper right
        local_wkt_extent << local_extent.MinX << " " << local_extent.MaxY << ", ";  // upper left
        local_wkt_extent << local_extent.MinX << " " << local_extent.MinY;         // lower left (again)
        local_wkt_extent << ")";
        */
        GEOSGeometry* geom_extent = OGREnvelope2GEOSGeom(local_extent);
        string local_wkt_extent = geos.as_string(geom_extent);
        GEOSGeom_destroy_r(geos.handle(), geom_extent);

#if 0
        #pragma omp critical
        {
            cout << "Thread " << omp_get_thread_num() << " " <<
                "MinX: " << local_extent.MinX << " "
                "MaxX: " << local_extent.MaxX << " "
                "MinY: " << local_extent.MinY << " "
                "MaxY: " << local_extent.MaxY <<
                endl;
            cout << "local_width: " << local_width << ", local_height: " << local_height << endl;
            cout << "Thread " << omp_get_thread_num() << " -- row: " << row << ", col: " << col << endl;
            cout << "Thread " << omp_get_thread_num() << " " << local_wkt_extent.str() << endl << endl;
        }
#endif

        zoneInfo* zone = new zoneInfo(local_extent, 0);

        // int count = 0;
        // OGREnvelope envelope;
        // for (int l1idx = 0; l1idx < lines.size(); ++l1idx) {
        //     lines[l1idx]->getEnvelope(&envelope);
        //     if (local_extent.Contains(envelope)) {
        //         ++count;
        //         zone->second.push_back(lines[l1idx]);
        //     }
        // }

        /**
         * This will select a line count of lines within our zone.
         */
        ostringstream oss;
        oss << "SELECT COUNT(1) FROM way WHERE "
            << "ST_SetSRID(ST_GeomFromText('" << local_wkt_extent << "'), 3395) ~ line2d_m";

        int nlines = 0;

        PGresult* res = db.query(oss.str().c_str());
        zone->second = atoi(PQgetvalue(res, 0, 0));

        if (zone->second > 0) {
            //cout << omp_get_thread_num() << ": " << oss.str() << endl;
            //cout << omp_get_thread_num() << ": got " << zone->second << " lines." << endl;
        }

        PQclear(res);

        sum_nb_lines += zone->second;

        #pragma omp critical
        zones.push_back(zone);
    }

    GEOSFree_r(hdl, hex_extent);

    if (maxdepth == 1) {
        return;
    }

    vector<int>                 to_del;
    vector< vector<zoneInfo*> > to_add;

    for (int zIdx = 0; zIdx < zones.size(); ++zIdx) {
        zoneInfo* zone = zones[zIdx];
        if (zone->second > 15000) {
            vector<zoneInfo*> subzones;
            GEOSGeometry* g = OGREnvelope2GEOSGeom(zone->first);
            prepare_zones(geos, g, subzones, maxdepth-1);
            GEOSGeom_destroy_r(hdl, g);

            to_del.insert(to_del.begin(), zIdx);
            to_add.push_back(subzones);
        }
    }

    for (int idx : to_del) {
        delete zones[idx];
        zones.erase(zones.begin() + idx);
    }
    for (auto subzones : to_add) {
        zones.insert(zones.end(), subzones.begin(), subzones.end());
    }
}

void write_zones(const std::string& filename, const vector<zoneInfo*>& zones, bool overwrite)
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

    OGRFieldDefn oField("color", OFTInteger);
    assert (z_layer->CreateField(&oField) == OGRERR_NONE);

    cout <<
        "MinX: " << zones[0]->first.MinX << " "
        "MaxX: " << zones[0]->first.MaxX << " "
        "MinY: " << zones[0]->first.MinY << " "
        "MaxY: " << zones[0]->first.MaxY <<
        endl;

    //return 0;

    for (int zIdx = 0; zIdx < zones.size(); ++zIdx) {
        OGREnvelope& env = zones[zIdx]->first;

        OGRFeature* poFeature = OGRFeature::CreateFeature(z_layer->GetLayerDefn());
        poFeature->SetField("color", rand()%0xFFFF);
        OGRPolygon* poly = new OGRPolygon();

        // per the specs, points must be added counter-clockwise
        OGRLinearRing* ring = new OGRLinearRing();
        ring->addPoint(zones[zIdx]->first.MinX, zones[zIdx]->first.MinY);
        ring->addPoint(zones[zIdx]->first.MaxX, zones[zIdx]->first.MinY);
        ring->addPoint(zones[zIdx]->first.MaxX, zones[zIdx]->first.MaxY);
        ring->addPoint(zones[zIdx]->first.MinX, zones[zIdx]->first.MaxY);
        ring->addPoint(zones[zIdx]->first.MinX, zones[zIdx]->first.MinY);    // this closes the ring

        poly->addRingDirectly(ring);

        // the feature takes ownership of the polygon, don't delete it.
        poFeature->SetGeometryDirectly(poly);

        assert (z_layer->CreateFeature(poFeature) == OGRERR_NONE);

        OGRFeature::DestroyFeature(poFeature);
    }

    GDALClose(zonesDS);
}

GEOSGeometry* world_geom()
{
    OGREnvelope env;

    #if 0
    env.MinX = -180.;
    env.MaxX = 180.;
    env.MinY = -90.;
    env.MaxY = 90.;
    #endif // 0

    env.MinX = -20026376.39;
    env.MaxX = 20026376.39;
    env.MinY = -15496570.74;
    env.MaxY = 18764656.23;

    return OGREnvelope2GEOSGeom(env);

    // GEOSCoordSequence* corners = GEOSCoordSeq_create_r(hdl, 5, 2);
    // assert (corners != NULL);
    //
    // GEOSCoordSeq_setX_r(hdl, corners, 0, -180.);
    // GEOSCoordSeq_setY_r(hdl, corners, 0, 90.);
    //
    // GEOSCoordSeq_setX_r(hdl, corners, 1, -180.);
    // GEOSCoordSeq_setY_r(hdl, corners, 1, -90.);
    //
    // GEOSCoordSeq_setX_r(hdl, corners, 2, 180.);
    // GEOSCoordSeq_setY_r(hdl, corners, 2, -90.);
    //
    // GEOSCoordSeq_setX_r(hdl, corners, 3, 180.);
    // GEOSCoordSeq_setY_r(hdl, corners, 3, 90.);
    //
    // GEOSCoordSeq_setX_r(hdl, corners, 4, -180.);
    // GEOSCoordSeq_setY_r(hdl, corners, 4, 90.);
    //
    // GEOSGeometry* shell = GEOSGeom_createLinearRing_r(hdl, corners);
    // assert (shell != NULL);
    //
    // GEOSGeometry* world_polygon = GEOSGeom_createPolygon_r(
    //     hdl,
    //     shell,
    //     NULL,
    //     0
    // );
    //
    // assert (world_polygon != NULL);
    //
    // GEOSSetSRID_r(hdl, world_polygon, 4326);
    // return world_polygon;
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

} // namespace cma
