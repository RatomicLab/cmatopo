#ifndef __CMA_PG_H
#define __CMA_PG_H

#include <string>
#include <libpq-fe.h>

namespace cma {

class PG
{
public:
    PG(const std::string& connect_str);
    ~PG();

    bool connected() const;
    PGresult* query(const std::string& sql, bool single_row_mode=false);
    PGresult* next_result();
    bool success(PGresult* res) const;

private:
    PGconn* _conn = NULL;
};

} // namespace cma

#endif // __CMA_PG_H
