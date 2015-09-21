#include <pg.h>

#include <iostream>

using namespace cma;
using namespace std;

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
    return _conn != NULL;
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
