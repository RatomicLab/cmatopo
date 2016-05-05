# cmatopo

A parallel, [PostGIS](http://postgis.net/) compatible, C++ topology creator.

## How it works

* The input is a PostgreSQL table/column (way.line2d_m) of the world's road network (linestrings, projected in meters -- SRID=3395) that needs to be added to a topology.
* Using MPI, it will process topology regions and will merge the result back into a single .ser file.
* This file can be converted to PostGIS compatible CSV files using the serdump utility.

## Dependencies

* PostgreSQL 9.4.5
* PostGIS 2.1.8, in particular the liblwgeom library and it's sources
* GEOS 3.6.0 (dev)
* GDAL 2.0.1
* Boost 1.60
* OpenMPI 1.8.8

# Compiling

This software was tested using the Intel C++ Compiler 14.0.

* Ajust path to liblwgeom (in-source) in the Makefile
* make -j8
