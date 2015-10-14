#!/usr/bin/env python

# cmatopo node output format:
#   1 | 2147483647 | POINT (-7819709.9884425168856978 6107814.9348522918298841)
#
# postgresql/postgis output format:
#   1,,POINT(-7819709.98844252 6107814.93485229)
#
#   postgis=# \f ','
#   Field separator is ",".
#   postgis=# \a
#   Output format is unaligned.
#   postgis=# \o 'postgis_output.csv'
#   postgis=# select node_id, containing_face, ST_AsText(geom) from way_topo.node order by node_id;

import sys
import itertools

### COMPARE NODES

cmatopo = sys.argv[1]
postgis = sys.argv[2]

with open(cmatopo, 'r') as f:
    cmatopo_lines = f.readlines()
with open(postgis, 'r') as f:
    postgis_lines = f.readlines()

if len(cmatopo_lines) != len(postgis_lines):
    print "ERR: different number of nodes. cmatopo: {0}, postgis: {1}.".format(len(cmatopo_lines), len(postgis_lines))

cmas = []
pgiss = []
for i in range(min(len(cmatopo_lines), len(postgis_lines))):
    cma = map(lambda e: e.strip(), cmatopo_lines[i].split('|'))
    cmapt = cma[2][7:-1].split(' ')
    cmas.append([int(cma[0]), map(lambda c: c[:c.find('.')+5], cmapt)])  # cma == [5476, ['-7737020.7574863', '6061074.4988690']]

    pgis = map(lambda e: e.strip(), postgis_lines[i].split(','))
    pgispt = pgis[2][6:-1].split(' ')
    pgiss.append([int(pgis[0]),  map(lambda c: c[:c.find('.')+5], pgispt)])

for i in range(len(cmas)):
    cma = cmas[i]
    pgis = pgiss[i]
    assert (cma[0] == pgis[0])

    if cma[1][0] != pgis[1][0] or cma[1][1] != pgis[1][1]:
        try:
            print "ERR: Point with id {0} doesn't match. ({1} != {2})".format(cma[0], "POINT({0} {1})".format(cma[1][0], cma[1][1]), "POINT({0} {1})".format(pgis[1][0], pgis[1][1]))
        except IOError:  # head, tail, less, ...
            break

### COMPARE EDGES

cmatopo = sys.argv[3]
postgis = sys.argv[4]

with open(cmatopo, 'r') as f:
    cmatopo_lines = f.readlines()
with open(postgis, 'r') as f:
    postgis_lines = f.readlines()

if len(cmatopo_lines) != len(postgis_lines):
    print "ERR: different number of edges. cmatopo: {0}, postgis: {1}.".format(len(cmatopo_lines), len(postgis_lines))

cmas = []
pgiss = []
for i in range(min(len(cmatopo_lines), len(postgis_lines))):
    cma = map(lambda e: e.strip(), cmatopo_lines[i].split('|'))
    cmaedge = list(itertools.chain(*map(lambda x: x.split(' '), map(lambda x: x.strip(), cma[9][12:-1].split(',')))))
    cmas.append([int(cma[0]), map(lambda c: c[:c.find('.')+5], cmaedge)])

    pgis = map(lambda e: e.strip(), postgis_lines[i].split('|'))
    pgisedge = list(itertools.chain(*map(lambda x: x.split(' '), map(lambda x: x.strip(), pgis[9][11:-1].split(',')))))
    pgiss.append([int(pgis[0]),  map(lambda c: c[:c.find('.')+5], pgisedge)])

for i in range(len(cmas)):
    cma = cmas[i]
    pgis = pgiss[i]
    assert (cma[0] == pgis[0])

    if len(cma[1]) != len(pgis[1]):
        print "ERR: Edge with id {0} doesn't have the same number of points ({1} != {2}).".format(cma[0], "LINESTRING({0})".format(" ".join(cma[1])), "LINESTRING({0})".format(" ".join(pgis[1])))

    for j in range(len(cma[1])):
        if cma[1][j] != pgis[1][j]:
            try:
                print "ERR: Edge with id {0} doesn't match. ({1} != {2})".format(cma[0], "LINESTRING({0})".format(" ".join(cma[1])), "LINESTRING({0})".format(" ".join(pgis[1])))
            except IOError:  # head, tail, less, ...
                break
