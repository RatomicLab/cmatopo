#include <zones.h>

#include <omp.h>
#include <assert.h>

#include <string>
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
const int cols = 4;

void prepare_zones(const linesV& lines, const OGREnvelope& extent, vector<zoneInfo*>& zones, int maxdepth)
{
    assert (maxdepth >= 1);
    assert (!omp_in_parallel());

    int nThreads = get_nb_threads();
    double width  = fabs(extent.MaxX - extent.MinX);
    double height = fabs(extent.MaxY - extent.MinY);

#if 0
    cout <<
        "MinX: " << extent.MinX << " "
        "MaxX: " << extent.MaxX << " "
        "MinY: " << extent.MinY << " "
        "MaxY: " << extent.MaxY <<
        endl;
#endif

    assert (rows*cols == nThreads);
    double local_width  = width / rows;
    double local_height = height / cols;

#if 0
    cout << "local_width: " << local_width << ", local_height: " << local_height << endl;
#endif

    double sum_nb_lines = 0.;

    #pragma omp parallel reduction(+:sum_nb_lines)
    {
        int row = int(floor(omp_get_thread_num() / rows));
        int col = omp_get_thread_num() % cols;

#if 0
        #pragma omp critical
        cout << "Thread " << omp_get_thread_num() << " -- row: " << row << ", col: " << col << endl;
#endif

        OGREnvelope local_extent;
        local_extent.MinX = extent.MinX + col * local_width;
        local_extent.MaxX = local_extent.MinX + local_width;
        local_extent.MinY = extent.MinY + row * local_height;
        local_extent.MaxY = local_extent.MinY + local_height;

#if 0
        #pragma omp critical
        cout <<
            "Local extent for " << omp_get_thread_num() << ": " <<
            " MinX: " << local_extent.MinX << " "
            " MaxX: " << local_extent.MaxX << " "
            " MinY: " << local_extent.MinY << " "
            " MaxY: " << local_extent.MaxY <<
            endl;
#endif

        zoneInfo* zone = new zoneInfo(local_extent, linesV());

        int count = 0;
        OGREnvelope envelope;
        for (int l1idx = 0; l1idx < lines.size(); ++l1idx) {
            lines[l1idx]->getEnvelope(&envelope);
            if (local_extent.Contains(envelope)) {
                ++count;
                zone->second.push_back(lines[l1idx]);
            }
        }

        sum_nb_lines += zone->second.size();

        #pragma omp critical
        zones.push_back(zone);

        //#pragma omp critical
        //cout << "Thread " << omp_get_thread_num() << " -- road count: " << zone->second.size() << endl;
    }

    double mean_zone_pop = sum_nb_lines / nThreads;

    if (mean_zone_pop < 4000 || maxdepth == 1) {
        return;
    }

#if 0
    cout << "Mean number of routes per zone: " << mean_zone_pop << endl;
#endif

    vector<int>                 to_del;
    vector< vector<zoneInfo*> > to_add;

    for (int zIdx = 0; zIdx < zones.size(); ++zIdx) {
        zoneInfo* zone = zones[zIdx];
        if (zone->second.size() > mean_zone_pop) {
            vector<zoneInfo*> subzones;
            prepare_zones(zone->second, zone->first, subzones, maxdepth-1);

            to_del.insert(to_del.begin(), zIdx);
            to_add.push_back(subzones);
        }
    }

    for (int idx : to_del) zones.erase(zones.begin() + idx);
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

    OGRSFDriver* shpDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
    assert (shpDriver != NULL);
    OGRDataSource *zonesDS = shpDriver->CreateDataSource("zones.shp", NULL);
    assert (zonesDS != NULL);
    OGRLayer* z_layer = zonesDS->CreateLayer("zones", NULL, wkbPolygon, NULL);
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

    OGRDataSource::DestroyDataSource(zonesDS);
}

} // namespace cma
