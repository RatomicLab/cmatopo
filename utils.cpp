#include <utils.h>

#include <omp.h>

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
    // TODO: do something with the message here eventually
    return;
}

} // namespace cma
