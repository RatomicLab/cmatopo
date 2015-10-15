#ifndef __CMA_UTILS_H
#define __CMA_UTILS_H

#include <geos_c.h>

#include <string>
#include <cassert>

namespace cma {

extern GEOSContextHandle_t hdl;

int get_nb_threads();
void geos_message_function(const char *fmt, ...);

class GEOSHelper
{
public:
    GEOSHelper() {
        hdl = initGEOS_r(geos_message_function, geos_message_function);
        if (hdl) {
            wkbr = GEOSWKBReader_create_r(hdl);
            wktr = GEOSWKTReader_create_r(hdl);
            wkbw = GEOSWKBWriter_create_r(hdl);
            wktw = GEOSWKTWriter_create_r(hdl);

            GEOSWKTWriter_setRoundingPrecision_r(hdl, wktw, 8);
        }
    }

    ~GEOSHelper() {
        if (wkbr) GEOSWKBReader_destroy_r(hdl, wkbr);
        if (wktr) GEOSWKTReader_destroy_r(hdl, wktr);
        if (wkbw) GEOSWKBWriter_destroy_r(hdl, wkbw);
        if (wktw) GEOSWKTWriter_destroy_r(hdl, wktw);

        finishGEOS_r(hdl);
    }

    GEOSContextHandle_t handle() const {
        return hdl;
    }

    GEOSWKBReader* reader() {
        return wkbr;
    }

    GEOSWKTReader* text_reader() {
        return wktr;
    }

    GEOSWKBWriter* writer() {
        return wkbw;
    }

    GEOSWKTWriter* text_writer() {
        return wktw;
    }

    std::string as_string(const GEOSGeometry* geom) {
        assert (geom);
        char* wkt_c = GEOSWKTWriter_write_r(hdl, text_writer(), geom);
        std::string wkt(wkt_c);
        GEOSFree_r(hdl, wkt_c);
        return wkt;
    }

private:
    GEOSWKBReader* wkbr = NULL;
    GEOSWKTReader* wktr = NULL;
    GEOSWKBWriter* wkbw = NULL;
    GEOSWKTWriter* wktw = NULL;
};

} // namespace cma

#endif // __CMA_UTILS_H
