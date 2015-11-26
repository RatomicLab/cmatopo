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
    GEOSGeometry* get_line(int id);

    bool get_lines(
        const GEOSGeometry* envelope,
        linesV& lines,
        bool within=true,
        int limit=-1);

    bool get_line_ids(
        const GEOSGeometry* envelope,
        std::set<int>& line_ids,
        bool within=true,
        int limit=-1);

    bool get_lines(
        const std::set<int> lineIds,
        linesV& lines);

private:
    std::string _build_query(
        const GEOSGeometry* geom,
        bool within,
        int limit,
        bool id);

    PGconn* _conn = NULL;
};

} // namespace cma

#endif // __CMA_PG_H
