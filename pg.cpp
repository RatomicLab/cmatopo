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

PG::PG(const std::string& connect_str)
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

PGresult* PG::query(const std::string& sql, bool single_row_mode)
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

bool PG::get_lines_within(OGREnvelope& envelope, linesV& lines)
{
    GEOSGeometry* geom = OGREnvelope2GEOSGeom(envelope);
    bool ret = get_lines_within(geom, lines);
    GEOSGeom_destroy_r(hdl, geom);
}

bool PG::get_lines_within(const GEOSGeometry* geom, linesV& lines)
{
    GEOSWKTWriter* wkt_writer = GEOSWKTWriter_create_r(hdl);

    char* hex = GEOSWKTWriter_write_r(hdl, wkt_writer, geom);
    ostringstream oss;
    oss << "SELECT line2d_m FROM way WHERE ST_SetSRID('" << hex << "'::geometry, 3395) ~ line2d_m order by id";
    cout << oss.str() << endl;
    GEOSFree_r(hdl, hex);
    GEOSWKTWriter_destroy_r(hdl, wkt_writer);

    PGresult* res = query(oss.str().c_str());

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

} // namespace cma
