#include <utils.h>

#include <omp.h>
#include <cstdio>
#include <stdarg.h>
#include <iostream>

using namespace std;

namespace cma {

int get_nb_threads()
{
    int nThreads = 1;
    #pragma omp parallel
    if (omp_get_thread_num() == 0) nThreads = omp_get_num_threads();
    return nThreads;
}

void geos_message_function(const char *fmt, ...)
{
    va_list ap, ap2;
    va_start(ap, fmt);
    printf(fmt, ap);
    va_end (ap);
}

void GEOSHelper::print_geom(const GEOSGeometry* geom)
{
    cout << as_string(geom) << endl;
}

} // namespace cma
