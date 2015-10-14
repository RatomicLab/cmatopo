#include <types.h>

#include <cassert>

namespace cma {
    extern GEOSContextHandle_t hdl;
}

namespace cma {

const GEOSPreparedGeometry* geom_container::prepared()
{
    assert (geom);
    if (!_prepared) {
        _prepared = GEOSPrepare_r(hdl, geom);
    }
    return _prepared;
}

const GEOSGeometry* geom_container::envelope()
{
    assert (geom);
    if (!_envelope) {
        _envelope = GEOSEnvelope_r(hdl, geom);
    }
    return _envelope;
}

geom_container::~geom_container()
{
    if (_prepared) {
        GEOSPreparedGeom_destroy_r(hdl, _prepared);
    }

    if (_envelope) {
        GEOSGeom_destroy_r(hdl, _envelope);
    }

    if (geom) {
        GEOSGeom_destroy_r(hdl, geom);
    }
}

bool geom_container::intersects(const GEOSGeometry* geom)
{
    return GEOSIntersects_r(hdl, envelope(), geom) == 1;
}

bool node::intersects(const GEOSGeometry* geom)
{
    if (GEOSGeomTypeId_r(hdl, geom) == GEOS_POINT) {
        double x1, y1, x2, y2;
        GEOSGeomGetX_r(hdl, this->geom, &x1);
        GEOSGeomGetY_r(hdl, this->geom, &y1);
        GEOSGeomGetX_r(hdl, this->geom, &x2);
        GEOSGeomGetY_r(hdl, this->geom, &y2);
        return x1 == x2 && y1 == y2;
    }
    return geom_container::intersects(geom);
}

} // namespace cma
