MPICXX ?= mpicxx

use_parallel_hdf5 = no
use_healpix = yes

PLUGINS = ./plugins
KEYVALUE_PATH = $(PLUGINS)/keyValue
HEALPIX_PATH = $(PLUGINS)/HEALPix/
QHULL_PATH = $(PLUGINS)/QHull/

HDF5_PLUGIN = $(PLUGINS)/hdf5

HDF5_PATH = /Users/jppaarde/sw/hdf5-1.6.10/installed
HDF5_FLAGS ?= -I$(HDF5_PATH)/include
HDF5_LIBS  ?= -L$(HDF5_PATH)/lib -lhdf5 -lz
HDF5_PLUGIN_FLAGS = -I$(HDF5_PLUGIN)

GSL_PATH = /usr/local/scisoft
GSL_FLAGS ?= -I$(GSL_PATH)/
GSL_LIBS ?= -L$(GSL_PATH)/lib -lgslcblas -lgsl 

SRC=./src/
OBJ=./obj/
LIB=./lib/
BIN=./bin/

CXXCFLAGS_NO_C = -Wall -ansi -O2 -g -ffast-math -fomit-frame-pointer
CXXCFLAGS := $(CXXCFLAGS_NO_C) -c
export CXXCFLAGS

MACRODEFS = -DH5_USE_16_API
MYLIBS = -L$(QHULL_PATH) -lqhull $(GSL_LIBS)  $(HDF5_LIBS)
MYINCLUDE = -I$(KEYVALUE_PATH) -I$(QHULL_PATH) $(HDF5_PLUGIN_FLAGS) $(HDF5_FLAGS) \
            -I$(PLUGINS)/octree -I$(PLUGINS)/hilbert -I$(SRC) \
            -I$(PLUGINS)/unit_sphere_tess/include $(GSL_FLAGS)
ifeq ($(use_parallel_hdf5),yes)
  MACRODEFS += -D HDF5_PARALLEL
  OBJS = $(OBJ)h5w_parallel.o
else
  OBJS = $(OBJ)h5w_serial.o
endif

ifeq ($(use_healpix),yes)
  MACRODEFS += -DHEAL_PIX
  MYINCLUDE += -I$(HEALPIX_PATH)/include
  MYLIBS += -L$(HEALPIX_PATH)/lib  -lhealpix_cxx
endif


OBJS += $(OBJ)SimpleX.o  $(OBJ)Structs.o $(OBJ)Common.o $(OBJ)rates.o	 \
        $(OBJ)keyValue.o $(OBJ)hilbert.o $(OBJ)tree_structures.o $(OBJ)Main.o


all: HEALPix QHull SimpleX

amuse: HEALPix QHull SimpleX amuse_interface

obj:
	mkdir obj

HEALPix:
	@cd $(HEALPIX_PATH)/src ; \
	make all

QHull:
	@cd  $(QHULL_PATH) ; \
	make all

SimpleX: HEALPix $(OBJS)
	$(MPICXX) $(OBJS) -o $(BIN)SimpleX $(CXXCFLAGS_NO_C) -lm $(MYLIBS)

$(OBJ)keyValue.o: $(KEYVALUE_PATH)/configfile.cpp $(KEYVALUE_PATH)/configfile.h
	$(MPICXX) -o$(OBJ)keyValue.o -c $(KEYVALUE_PATH)/configfile.cpp $(CXXCFLAGS) 

$(OBJ)tree_structures.o: $(PLUGINS)/octree/tree_structures.cpp $(PLUGINS)/octree/tree_structures.h $(SRC)Structs.h
	$(MPICXX) -o$(OBJ)tree_structures.o $(MACRODEFS) -c $(PLUGINS)/octree/tree_structures.cpp $(CXXCFLAGS) $(MYINCLUDE)

$(OBJ)hilbert.o: $(PLUGINS)/hilbert/hilbert.c $(PLUGINS)/hilbert/hilbert.h
	$(MPICXX) -o$(OBJ)hilbert.o -c $(PLUGINS)/hilbert/hilbert.c $(CXXCFLAGS) $(MYINCLUDE)

$(OBJ)Structs.o: $(SRC)Structs.cpp $(SRC)Structs.h
	$(MPICXX) -o$(OBJ)Structs.o -c $(SRC)Structs.cpp $(CXXCFLAGS)

$(OBJ)Common.o: $(SRC)Common.cpp $(SRC)Common.h
	$(MPICXX) -o$(OBJ)Common.o -c $(SRC)Common.cpp $(CXXCFLAGS)

$(OBJ)rates.o: $(SRC)rates.h $(SRC)rates.cpp
	$(MPICXX) -o$(OBJ)rates.o -c $(SRC)rates.cpp $(CXXCFLAGS)

$(OBJ)SimpleX.o: $(SRC)SimpleX.cpp $(SRC)SimpleX.h $(SRC)Structs.h $(SRC)Common.h $(SRC)rates.h
	$(MPICXX) -o$(OBJ)SimpleX.o $(MYINCLUDE) $(MACRODEFS) -c $(SRC)SimpleX.cpp $(CXXCFLAGS) 

$(OBJ)Main.o: $(SRC)Main.cpp $(SRC)Main.h $(SRC)SimpleX.h $(SRC)Common.h
	$(MPICXX) -o$(OBJ)Main.o $(MYINCLUDE) -c $(SRC)Main.cpp $(CXXCFLAGS) $(MACRODEFS)

$(OBJ)h5w_parallel.o: $(HDF5_PLUGIN)/h5w_parallel.cpp $(HDF5_PLUGIN)/h5w_parallel.h $(HDF5_PLUGIN)/array_class.h
	$(MPICXX) -o$(OBJ)h5w_parallel.o -c $(HDF5_PLUGIN)/h5w_parallel.cpp $(CXXCFLAGS) $(MACRODEFS)

$(OBJ)h5w_serial.o: $(HDF5_PLUGIN)/h5w_serial.cpp $(HDF5_PLUGIN)/h5w_serial.h $(HDF5_PLUGIN)/array_class.h
	$(MPICXX) -o$(OBJ)h5w_serial.o -c $(HDF5_PLUGIN)/h5w_serial.cpp $(CXXCFLAGS) $(MACRODEFS) $(HDF5_FLAGS) $(HDF5_PLUGIN_FLAGS)

amuse_interface: $(OBJ)amuse_interface.o 
amuse_interface_gl: $(OBJ)amuse_interface_gl.o 

$(OBJ)amuse_interface.o: $(SRC)amuse_interface.cpp $(SRC)amuse_interface.h $(SRC)SimpleX.cpp $(SRC)SimpleX.h $(SRC)Structs.h $(SRC)Common.h $(PLUGINS)/octree/tree_structures.cpp
	$(MPICXX) -o$(OBJ)amuse_interface.o $(MYINCLUDE) $(MACRODEFS) -c $(SRC)amuse_interface.cpp $(CXXCFLAGS)

$(OBJ)amuse_interface_gl.o: $(SRC)amuse_interface_gl.cpp $(OBJ)amuse_interface.o $(SRC)SimpleX.cpp $(SRC)SimpleX.h $(SRC)Structs.h $(SRC)Common.h $(PLUGINS)/octree/tree_structures.cpp
	$(MPICXX) -o$(OBJ)amuse_interface_gl.o $(MYINCLUDE) $(MACRODEFS) -c $(SRC)amuse_interface_gl.cpp $(CXXCFLAGS)

clean:
	rm -rf $(OBJ)*.o $(BIN)SimpleX $(BIN)SimpleX.exe

veryclean: clean
	@cd $(HEALPIX_PATH)/src ; \
	make clean
	@cd  $(QHULL_PATH) ; \
	make clean

distclean: clean veryclean
