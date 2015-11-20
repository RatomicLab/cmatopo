#include <pg.h>

#include <cassert>
#include <sstream>
#include <iostream>

#include <zones.h>

using namespace cma;
using namespace std;

namespace cma {
    extern GEOSContextHandle_t hdl;
}

namespace cma {

PG::PG(const string& connect_str)
{
    _conn = PQconnectdb(connect_str.c_str());
}

PG::~PG()
{
    if (_conn) {
        PQfinish(_conn);
    }
}

bool PG::connected() const
{
    return _conn != NULL && PQstatus(_conn) == CONNECTION_OK;
}

PGresult* PG::query(const string& sql, bool single_row_mode)
{
    if (!single_row_mode) {
        return PQexec(_conn, sql.c_str());
    }

    int qid = PQsendQuery(_conn, sql.c_str());
    if (PQsetSingleRowMode(_conn) != 1) {
        cerr << "Warning: could not set single-row mode for query: " << sql << endl;
    }
    return next_result();
}

PGresult* PG::next_result()
{
    return PQgetResult(_conn);
}

bool PG::success(PGresult* res) const
{
    if (!res) {
        return false;
    }

    int status = PQresultStatus(res);
    return (status == PGRES_TUPLES_OK || status == PGRES_COMMAND_OK);
}

int PG::get_line_count()
{
    ostringstream oss;
    oss << "SELECT count(1) FROM way;";

    PGresult* res = query(oss.str().c_str());
    if (!success(res)) {
        PQclear(res);
        return -1;
    }

    char* count = PQgetvalue(res, 0, 0);
    PQclear(res);
    return atoi(count);
}

bool PG::get_lines(
    const GEOSGeometry* geom,
    linesV& lines,
    bool within,
    int limit)
{
    string q = _build_query(geom, within, limit, false);

    PGresult* res = query(q.c_str());

    if (!success(res)) {
        PQclear(res);
        return false;
    }

    char* line2d;
    GEOSWKBReader* wkb_reader = GEOSWKBReader_create_r(hdl);
    for (int i = 0; i < PQntuples(res); i++) {
        line2d = PQgetvalue(res, i, 0);
        GEOSGeometry* line = GEOSWKBReader_readHEX_r(hdl, wkb_reader, (const unsigned char*)line2d, PQgetlength(res, i, 0));
        assert (line != NULL);
        lines.push_back(line);
    }
    GEOSWKBReader_destroy_r(hdl, wkb_reader);

    PQclear(res);

    return true;
}

bool PG::get_line_ids(
    const GEOSGeometry* envelope,
    set<int>& line_ids,
    bool within,
    int limit)
{
    string q = _build_query(envelope, within, limit, true);

    PGresult* res = query(q.c_str());

    if (!success(res)) {
        PQclear(res);
        return false;
    }

    char* id;
    for (int i = 0; i < PQntuples(res); i++) {
        id = PQgetvalue(res, i, 0);
        assert (id != NULL);
        line_ids.insert(atoi(id));
    }

    PQclear(res);

    return true;
}

string PG::_build_query(
    const GEOSGeometry* geom,
    bool within,
    int limit,
    bool id)
{
    GEOSWKTWriter* wkt_writer = GEOSWKTWriter_create_r(hdl);
    char* hex = GEOSWKTWriter_write_r(hdl, wkt_writer, geom);

    ostringstream oss;
    oss << "SELECT " << (id ? "id" : "line2d_m") << " FROM way WHERE ";
    if (!within) oss << " NOT ";
    oss << "ST_SetSRID('" << hex << "'::geometry, 3395) ~ line2d_m ORDER BY id";
    if (limit > 0) oss << " LIMIT " << limit;

    GEOSFree_r(hdl, hex);
    GEOSWKTWriter_destroy_r(hdl, wkt_writer);

    return oss.str();
}

} // namespace cma
