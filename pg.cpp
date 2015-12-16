#include <pg.h>

#include <string>
#include <cassert>
#include <sstream>
#include <iostream>

#include <boost/algorithm/string/join.hpp>

#include <zones.h>

using namespace cma;
using namespace std;

using namespace boost::algorithm;

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

    int count = atoi(PQgetvalue(res, 0, 0));
    PQclear(res);
    return count;
}

GEOSGeometry* PG::get_line(int id)
{
    ostringstream oss;
    oss << "SELECT line2d_m FROM way WHERE id=" << id;

    PGresult* res = query(oss.str().c_str());

    if (!success(res)) {
        PQclear(res);
        return nullptr;
    }

    char* line2d;
    GEOSWKBReader* wkb_reader = GEOSWKBReader_create_r(hdl);
    line2d = PQgetvalue(res, 0, 0);
    GEOSGeometry* line = GEOSWKBReader_readHEX_r(hdl, wkb_reader, (const unsigned char*)line2d, PQgetlength(res, 0, 0));
    assert (line != NULL);
    GEOSWKBReader_destroy_r(hdl, wkb_reader);

    return line;
}

bool PG::get_lines(
    const GEOSGeometry* geom,
    linesV& lines,
    bool within,
    int limit)
{
    string q = _build_query(geom, within, limit);

    PGresult* res = query(q.c_str());

    if (!success(res)) {
        PQclear(res);
        return false;
    }

    char* id;
    char* line2d;
    GEOSWKBReader* wkb_reader = GEOSWKBReader_create_r(hdl);
    for (int i = 0; i < PQntuples(res); i++) {
        id = PQgetvalue(res, i, 0);
        line2d = PQgetvalue(res, i, 1);
        GEOSGeometry* line = GEOSWKBReader_readHEX_r(hdl, wkb_reader, (const unsigned char*)line2d, PQgetlength(res, i, 1));
        if (line != NULL) {
            // assert (line != NULL);
            lines.push_back(make_pair(atoi(id), line));
        }
        else {
            cerr << "Line id " << id << " could not be converted to a GEOSGeometry..." << endl;
        }
    }
    GEOSWKBReader_destroy_r(hdl, wkb_reader);

    PQclear(res);

    return true;
}

bool PG::get_common_lines(
    const OGREnvelope& env1,
    const OGREnvelope& env2,
    linesV& lines,
    int limit)
{
    GEOSGeometry* g1 = OGREnvelope2GEOSGeom(env1);
    GEOSGeometry* g2 = OGREnvelope2GEOSGeom(env2);

    string env1_str = build_pg_geom(g1);
    string env2_str = build_pg_geom(g2);

    OGREnvelope m = env1;
    m.Merge(env2);

    GEOSGeometry* g3 = OGREnvelope2GEOSGeom(m);
    string m_str = build_pg_geom(g3);

    ostringstream oss;
    oss << "SELECT id, line2d_m FROM way WHERE "
        << m_str << " && line2d_m AND ST_Contains(" << m_str << ", line2d_m) AND "
        << env1_str << " && line2d_m AND NOT ST_Contains(" << env1_str << ", line2d_m) AND "
        << env2_str << " && line2d_m AND NOT ST_Contains(" << env2_str << ", line2d_m)";

    PGresult* res = query(oss.str().c_str());

    if (!success(res)) {
        PQclear(res);
        return false;
    }

    char* id;
    char* line2d;
    GEOSWKBReader* wkb_reader = GEOSWKBReader_create_r(hdl);
    for (int i = 0; i < PQntuples(res); i++) {
        id = PQgetvalue(res, i, 0);
        line2d = PQgetvalue(res, i, 1);
        GEOSGeometry* line = GEOSWKBReader_readHEX_r(hdl, wkb_reader, (const unsigned char*)line2d, PQgetlength(res, i, 1));
        assert (line != NULL);
        lines.push_back(make_pair(atoi(id), line));
    }
    GEOSWKBReader_destroy_r(hdl, wkb_reader);

    GEOSGeom_destroy_r(hdl, g1);
    GEOSGeom_destroy_r(hdl, g2);
    GEOSGeom_destroy_r(hdl, g3);
}

bool PG::get_line_ids(
    const GEOSGeometry* envelope,
    set<int>& line_ids,
    bool within,
    int limit)
{
    string q = _build_query(envelope, within, limit);

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

bool PG::get_lines(const std::set<int> lineIds, linesV& lines)
{
    vector<string> ids;
    transform(lineIds.begin(), lineIds.end(), back_inserter(ids), [](int a) {
        return to_string(a);
    });

    ostringstream oss;
    oss << "SELECT id, line2d_m FROM way WHERE id IN (" << join(ids, ",") << ")";

    PGresult* res = query(oss.str().c_str());

    if (!success(res)) {
        PQclear(res);
        return false;
    }

    char* id;
    char* line2d;
    GEOSWKBReader* wkb_reader = GEOSWKBReader_create_r(hdl);
    for (int i = 0; i < PQntuples(res); i++) {
        id = PQgetvalue(res, i, 0);
        line2d = PQgetvalue(res, i, 1);
        GEOSGeometry* line = GEOSWKBReader_readHEX_r(hdl, wkb_reader, (const unsigned char*)line2d, PQgetlength(res, i, 1));
        assert (line != NULL);
        lines.push_back(make_pair(atoi(id), line));
    }
    GEOSWKBReader_destroy_r(hdl, wkb_reader);

    PQclear(res);
    return true;
}

string PG::_build_query(
    const GEOSGeometry* geom,
    bool within,
    int limit)
{
    string geom_str = build_pg_geom(geom);

    ostringstream oss;
    oss << "SELECT id, line2d_m FROM way WHERE ";
    if (!within) oss << " NOT ";
    oss << " ST_Contains(" << geom_str << ", line2d_m) ";
    if (!within) oss << " AND ST_Intersects(" << geom_str << ", line2d_m) ";
    oss << " ORDER BY id";
    if (limit > 0) oss << " LIMIT " << limit;

    cout << oss.str() << endl;
    return oss.str();
}

string PG::build_pg_geom(const GEOSGeometry* geom) const
{
    GEOSWKTWriter* wkt_writer = GEOSWKTWriter_create_r(hdl);
    char* hex = GEOSWKTWriter_write_r(hdl, wkt_writer, geom);

    ostringstream oss_geom;
    oss_geom << " ST_SetSRID('" << hex << "'::geometry, 3395) ";

    GEOSFree_r(hdl, hex);
    GEOSWKTWriter_destroy_r(hdl, wkt_writer);

    return oss_geom.str();
}

} // namespace cma
