#include <fstream>
#include <iostream>
#include <boost/serialization/vector.hpp>

#include <geos_c.h>

#include <topology.h>

using namespace cma;
using namespace std;

namespace cma {
    GEOSContextHandle_t hdl;
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <filename.ser>" << endl;
        return 1;
    }

    ifstream ifs(argv[1]);
    if (!ifs.is_open()) {
        cerr << "Could not open file " << argv[1] << endl;
        return 1;
    }

    initGEOS(geos_message_function, geos_message_function);

    unique_ptr<GEOSHelper> geos(new GEOSHelper());

    Topology* t = new Topology(geos.get());
    boost::archive::binary_iarchive ia(ifs);
    ia >> *t;

    t->pg_output();

    cout << "Topology has been saved for PostgreSQL import." << endl
         << endl
         << "1. Please issue the following command to create your topology:" << endl
         << "\tselect topology.CreateTopology('way_topo', 3395, 1);" << endl
         << "\tselect topology.AddTopoGeometryColumn('way_topo', 'public', 'way', 'topo_geom', 'LINE');" << endl
         << endl
         << "2. Copy the node.csv, edge_data.csv, face.csv and relation.csv files over to your database server." << endl
         << endl
         << "3. Issue the following to import the generated topology (in order):" << endl
         << "\tCOPY way_topo.node FROM '/path/to/node.csv' WITH CSV;" << endl
         << "\tCOPY way_topo.face FROM '/path/to/face.csv' WITH CSV DELIMITER '|';" << endl
         << "\tCOPY way_topo.edge_data FROM '/path/to/edge_data.csv' WITH CSV DELIMITER '|';" << endl
         << "\tCOPY way_topo.relation FROM '/path/to/relation.csv' WITH CSV;" << endl;
         << endl
         << "4. Insert the content of the GetTopoGeom.sql file into your database (if not already done)." << endl
         << endl
         << "5. Issue the following to update the way.topo_geom column:" << endl
         << "\tcat /path/to/topo_geom.sql | psql <your_database>" << endl;

    delete t;

    finishGEOS();
    return 0;
}
