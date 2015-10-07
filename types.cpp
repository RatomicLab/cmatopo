#include <types.h>

namespace cma {
    extern GEOSContextHandle_t hdl;
}

namespace cma {

edge::~edge()
{
    if (geom) GEOSGeom_destroy_r(hdl, geom);
    if (_envelope) GEOSGeom_destroy_r(hdl, _envelope);
}

bool edge::intersects(const GEOSGeometry* geom)
{
    return GEOSIntersects_r(hdl, envelope(), geom) == 1;
}

const GEOSPreparedGeometry* edge::prepared()
{
    if (!_prepared) {
        _prepared = GEOSPrepare_r(hdl, geom);
    }
    return _prepared;
}

const GEOSGeometry* edge::envelope()
{
    if (!_envelope) {
        _envelope = GEOSEnvelope_r(hdl, geom);
    }
    return _envelope;
}

node::~node()
{
    if (geom) GEOSGeom_destroy_r(hdl, geom);
    if (_envelope) GEOSGeom_destroy_r(hdl, _envelope);
}

bool node::intersects(const GEOSGeometry* geom)
{
    return GEOSIntersects_r(hdl, envelope(), geom) == 1;
}

const GEOSPreparedGeometry* node::prepared()
{
    if (!_prepared) {
        _prepared = GEOSPrepare_r(hdl, geom);
    }
    return _prepared;
}

const GEOSGeometry* node::envelope()
{
    if (!_envelope) {
        _envelope = GEOSEnvelope_r(hdl, geom);
    }
    return _envelope;
}

face::~face()
{
    if (_envelope) GEOSGeom_destroy_r(hdl, _envelope);
}

bool face::intersects(const GEOSGeometry* geom)
{
    return GEOSIntersects_r(hdl, envelope(), geom) == 1;
}

const GEOSGeometry* face::envelope()
{
    if (!_envelope) {
        _envelope = GEOSEnvelope_r(hdl, mbr);
    }
    return _envelope;
}

} // namespace cma
