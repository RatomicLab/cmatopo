CXX := icpc

program_NAME := cmatopo
program_SRCS := $(filter-out serdump.cpp main.cpp, $(wildcard *.cpp))
program_OBJS := ${program_SRCS:.cpp=.o}

program_INCLUDE_DIRS := .
#program_INCLUDE_DIRS += /Users/laurent/software/elcapitan/postgis/2.1.8_intel15/include
program_INCLUDE_DIRS += /home/laurent/src/postgis-2.1.8/liblwgeom
#program_INCLUDE_DIRS += /Users/laurent/software/elcapitan/postgis/2.1.8_intel15_dbg/include
program_INCLUDE_DIRS += /usr/local/include

#program_LIBRARY_DIRS := /Users/laurent/software/elcapitan/postgis/2.1.8_intel15/lib
program_LIBRARY_DIRS := /home/laurent/src/postgis-2.1.8/liblwgeom
#program_LIBRARY_DIRS := /Users/laurent/software/elcapitan/postgis/2.1.8_intel15_dbg/lib
program_LIBRARY_DIRS += /usr/local/lib


program_LIBRARIES := mpfr boost_serialization-mt boost_filesystem-mt boost_system-mt boost_mpi-mt mpi lwgeom geos_c

CPPFLAGS =  -g -openmp -O3 -std=c++11 `gdal-config --cflags` `geos-config --cflags`
#CPPFLAGS =  -g -openmp -O0 -std=c++11 `gdal-config --cflags` `geos-config --cflags`
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))

LDFLAGS += `gdal-config --libs` `geos-config --libs`
LDFLAGS += -lpq
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))

.PHONY: all clean distclean

all: $(program_NAME) serdump

$(program_NAME): main.cpp $(program_OBJS)
	$(LINK.cc) $(program_OBJS) main.cpp -o $(program_NAME)

serdump: serdump.cpp $(program_OBJS)
	$(LINK.cc) $(program_OBJS) serdump.cpp -o serdump

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)
	@- $(RM) serdump

distclean: clean
