#ifndef __CMA_UTILS_H
#define __CMA_UTILS_H

#include <geos_c.h>

namespace cma {

int get_nb_threads();
void geos_message_function(const char *fmt, ...); 

} // namespace cma

#endif // __CMA_UTILS_H
