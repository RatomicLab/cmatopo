#ifndef __CMA_PG_H
#define __CMA_PG_H

#include <string>
#include <libpq-fe.h>
#include <ogrsf_frmts.h>

#include <types.h>

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
    int get_line_count();

    /**
     * Below are some geometry helpers.
     */
    bool get_lines_within(const GEOSGeometry* envelope, linesV& lines);
    bool get_lines_within(OGREnvelope& envelope, linesV& lines);

private:
    PGconn* _conn = NULL;
};

} // namespace cma

#endif // __CMA_PG_H
