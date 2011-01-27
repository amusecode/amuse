MPICXX ?= mpicxx

use_parallel_hdf5 = no
use_healpix = yes

PLUGINS = plugins
KEYVALUE_PATH = $(PLUGINS)/keyValue
HDF5_PATH = $(PLUGINS)/hdf5
HEALPIX_PATH = $(PLUGINS)/HEALPix/
QHULL_PATH = $(PLUGINS)/QHull/

SRC=./src/
OBJ=./obj/
LIB=./lib/
BIN=./bin/
CXXWARNFLAGS= -W -Wall -Wshadow -Wwrite-strings -Wredundant-decls -Woverloaded-virtual -Wcast-qual -Wcast-align -Wpointer-arith -Wconversion -Wold-style-cast -Wno-unknown-pragmas
CXXCFLAGS= $(CXXWARNFLAGS) -ansi -I$(INCDIR) -O2 -g0 -ffast-math -fomit-frame-pointer -c
#CXXCFLAGS = -O -g -Wall -W -Wno-deprecated
export CXXCFLAGS

MACRODEFS = -DH5_USE_16_API
MYLIBS = -L${QHULL_PATH} -lqhull -lgslcblas -lgsl -lhdf5 -lz
MYINCLUDE = -I${KEYVALUE_PATH} -I${QHULL_PATH} -I${HDF5_PATH}  \
            -I${PLUGINS}/octree -I${PLUGINS}/hilbert -I${SRC} \
            -I${PLUGINS}/unit_sphere_tess/include

ifeq ($(use_parallel_hdf5),yes)
  MACRODEFS += -D HDF5_PARALLEL
  OBJS = ${OBJ}h5w_parallel.o
else
  OBJS = ${OBJ}h5w_serial.o
endif

ifeq ($(use_healpix),yes)
  MACRODEFS += -DHEAL_PIX
  MYINCLUDE += -I${HEALPIX_PATH}/include
  MYLIBS += -L${HEALPIX_PATH}/lib  -lhealpix_cxx
endif


OBJS += ${OBJ}SimpleX.o  ${OBJ}Structs.o ${OBJ}Common.o \
        ${OBJ}keyValue.o ${OBJ}hilbert.o ${OBJ}tree_structures.o ${OBJ}Main.o


all: HEALPix QHull SimpleX

obj:
	mkdir obj
	
HEALPix:
	@cd $(HEALPIX_PATH)/src ; \
	make all

QHull:
	@cd  $(QHULL_PATH) ; \
	make all

SimpleX: HEALPix ${OBJS}
	${MPICXX} ${OBJS} -o ${BIN}SimpleX ${CXXCFLAGS} -lm $(MYLIBS)

${OBJ}keyValue.o: $(KEYVALUE_PATH)/configfile.cpp $(KEYVALUE_PATH)/configfile.h
	$(MPICXX) -o${OBJ}keyValue.o -c $(KEYVALUE_PATH)/configfile.cpp ${CXXCFLAGS} 

${OBJ}tree_structures.o: ${PLUGINS}/octree/tree_structures.cpp ${PLUGINS}/octree/tree_structures.h ${SRC}Structs.h
	${MPICXX} -o${OBJ}tree_structures.o ${MACRODEFS} -c ${PLUGINS}/octree/tree_structures.cpp ${CXXCFLAGS} ${MYINCLUDE}

${OBJ}hilbert.o: ${PLUGINS}/hilbert/hilbert.c ${PLUGINS}/hilbert/hilbert.h
	${MPICXX} -o${OBJ}hilbert.o -c ${PLUGINS}/hilbert/hilbert.c ${CXXCFLAGS} ${MYINCLUDE}

${OBJ}Structs.o: ${SRC}Structs.cpp ${SRC}Structs.h
	${MPICXX} -o${OBJ}Structs.o -c ${SRC}Structs.cpp ${CXXCFLAGS}

${OBJ}Common.o: ${SRC}Common.cpp ${SRC}Common.h
	${MPICXX} -o${OBJ}Common.o -c ${SRC}Common.cpp ${CXXCFLAGS}

${OBJ}SimpleX.o: ${SRC}SimpleX.cpp ${SRC}SimpleX.h ${SRC}Structs.h ${SRC}Common.h 
	${MPICXX} -o${OBJ}SimpleX.o ${MYINCLUDE} ${MACRODEFS} -c ${SRC}SimpleX.cpp ${CXXCFLAGS}

${OBJ}Main.o: ${SRC}Main.cpp ${SRC}Main.h ${SRC}SimpleX.h ${SRC}Common.h
	${MPICXX} -o${OBJ}Main.o ${MYINCLUDE} -c ${SRC}Main.cpp ${CXXCFLAGS} ${MACRODEFS}

${OBJ}h5w_parallel.o: ${HDF5_PATH}/h5w_parallel.cpp ${HDF5_PATH}/h5w_parallel.h ${HDF5_PATH}/array_class.h
	${MPICXX} -o${OBJ}h5w_parallel.o -c ${HDF5_PATH}/h5w_parallel.cpp ${CXXCFLAGS} ${MACRODEFS}

${OBJ}h5w_serial.o: ${HDF5_PATH}/h5w_serial.cpp ${HDF5_PATH}/h5w_serial.h ${HDF5_PATH}/array_class.h
	${MPICXX} -o${OBJ}h5w_serial.o -c ${HDF5_PATH}/h5w_serial.cpp ${CXXCFLAGS} ${MACRODEFS}

amuse_interface: ${OBJ}amuse_interface.o 
amuse_interface_gl: ${OBJ}amuse_interface_gl.o 

${OBJ}amuse_interface.o: ${SRC}amuse_interface.cpp ${SRC}SimpleX.cpp ${SRC}SimpleX.h ${SRC}Structs.h ${SRC}Common.h ${PLUGINS}/octree/tree_structures.cpp
	${MPICXX} -o${OBJ}amuse_interface.o ${MYINCLUDE} ${MACRODEFS} -c ${SRC}amuse_interface.cpp ${CXXCFLAGS}

${OBJ}amuse_interface_gl.o: ${SRC}amuse_interface_gl.cpp ${OBJ}amuse_interface.o ${SRC}SimpleX.cpp ${SRC}SimpleX.h ${SRC}Structs.h ${SRC}Common.h ${PLUGINS}/octree/tree_structures.cpp
	${MPICXX} -o${OBJ}amuse_interface_gl.o ${MYINCLUDE} ${MACRODEFS} -c ${SRC}amuse_interface_gl.cpp ${CXXCFLAGS}

clean:
	rm -rf ${OBJ}*.o ${BIN}SimpleX
	@cd $(HEALPIX_PATH)/src ; \
	make clean
	@cd  $(QHULL_PATH) ; \
	make clean
