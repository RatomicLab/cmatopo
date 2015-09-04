CXX := icpc

program_NAME := cmatopo
program_SRCS := $(wildcard *.cpp)
program_OBJS := ${program_SRCS:.cpp=.o}

program_INCLUDE_DIRS := . /usr/local/include /Users/laurent/Projets/CycleMapApp/postgis-2.1.8/liblwgeom
program_LIBRARY_DIRS := /Users/laurent/Projets/CycleMapApp/postgis-2.1.8/liblwgeom/.libs
program_LIBRARIES := mpfr boost_serialization-mt boost_filesystem-mt boost_system-mt boost_mpi-mt lwgeom geos_c

CPPFLAGS = -openmp -O3 -std=c++11 `gdal-config --cflags`
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))

LDFLAGS = `gdal-config --libs`
LDFLAGS += -lpq
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))

.PHONY: all clean distclean

all: $(program_NAME)

$(program_NAME): $(program_OBJS)
	$(LINK.cc) $(program_OBJS) -o $(program_NAME)

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)

distclean: clean

run:
	DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2015.1.108/compiler/lib:/Users/laurent/Projets/CycleMapApp/postgis-2.1.8/liblwgeom/.libs mpirun -np 4 ./test
