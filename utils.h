#ifndef __CMA_UTILS_H
#define __CMA_UTILS_H

#include <geos_c.h>

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
        }
    }
    
    ~GEOSHelper() {
        if (wkbr) GEOSWKBReader_destroy_r(hdl, wkbr);
        if (wktr) GEOSWKTReader_destroy_r(hdl, wktr);
        if (wkbw) GEOSWKBWriter_destroy_r(hdl, wkbw);

        finishGEOS_r(hdl);
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

private:
    GEOSWKBReader* wkbr = NULL;
    GEOSWKTReader* wktr = NULL;
    GEOSWKBWriter* wkbw = NULL;
};

} // namespace cma

#endif // __CMA_UTILS_H
