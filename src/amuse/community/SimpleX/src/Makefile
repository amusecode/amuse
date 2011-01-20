
use_parallel_hdf5 = no
use_healpix = yes

MPICXX = mpicxx
PLUGINS = plugins

KEYVALUE_PATH = $(PLUGINS)/keyValue
HDF5_PATH = $(PLUGINS)/hdf5
HEALPIX_PATH = $(PLUGINS)/HEALPix/Healpix_2.11/src/cxx/generic_gcc/

HDF5_INCL = -DH5_USE_16_API
HDF5_LIB =  -lhdf5 -lz 

SRC=./src/
OBJ=./obj/
LIB=./lib/
BIN=./bin/
FLAGS = -O -g -Wall -W -Wno-deprecated

MACRODEFS = 
MYLIBS = -L${PLUGINS}/QHull -lqhull -lgslcblas -lgsl
MYINCLUDE = -I${PLUGINS}/keyValue -I${PLUGINS}/QHull/src/ -I${PLUGINS}/hdf5  \
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


all: SimpleX

SimpleX: ${OBJS}
	${MPICXX} ${OBJS} -o ${BIN}SimpleX ${FLAGS} -lm $(MYLIBS) ${HDF5_LIB}

${OBJ}keyValue.o: $(KEYVALUE_PATH)/configfile.cpp $(KEYVALUE_PATH)/configfile.h
	$(MPICXX) -o${OBJ}keyValue.o -c $(KEYVALUE_PATH)/configfile.cpp ${FLAGS} 

${OBJ}tree_structures.o: ${PLUGINS}/octree/tree_structures.cpp ${PLUGINS}/octree/tree_structures.h ${SRC}Structs.h
	${MPICXX} -o${OBJ}tree_structures.o ${MACRODEFS} -c ${PLUGINS}/octree/tree_structures.cpp ${FLAGS} ${MYINCLUDE}

${OBJ}hilbert.o: ${PLUGINS}/hilbert/hilbert.c ${PLUGINS}/hilbert/hilbert.h
	${MPICXX} -o${OBJ}hilbert.o -c ${PLUGINS}/hilbert/hilbert.c ${FLAGS} ${MYINCLUDE}

${OBJ}Structs.o: ${SRC}Structs.cpp ${SRC}Structs.h
	${MPICXX} -o${OBJ}Structs.o -c ${SRC}Structs.cpp ${FLAGS}

${OBJ}Common.o: ${SRC}Common.cpp ${SRC}Common.h
	${MPICXX} -o${OBJ}Common.o -c ${SRC}Common.cpp ${FLAGS}

${OBJ}SimpleX.o: ${SRC}SimpleX.cpp ${SRC}SimpleX.h ${SRC}Structs.h ${SRC}Common.h 
	${MPICXX} -o${OBJ}SimpleX.o ${MYINCLUDE} ${MACRODEFS} -c ${SRC}SimpleX.cpp ${FLAGS} ${HDF5_INCL} 

${OBJ}Main.o: ${SRC}Main.cpp ${SRC}Main.h ${SRC}SimpleX.h ${SRC}Common.h
	${MPICXX} -o${OBJ}Main.o ${MYINCLUDE} -c ${SRC}Main.cpp ${FLAGS} ${HDF5_INCL} ${MACRODEFS}

${OBJ}h5w_parallel.o: ${HDF5_PATH}/h5w_parallel.cpp ${HDF5_PATH}/h5w_parallel.h ${HDF5_PATH}/array_class.h
	${MPICXX} -o${OBJ}h5w_parallel.o -c ${HDF5_PATH}/h5w_parallel.cpp ${FLAGS} ${HDF5_INCL}

${OBJ}h5w_serial.o: ${HDF5_PATH}/h5w_serial.cpp ${HDF5_PATH}/h5w_serial.h ${HDF5_PATH}/array_class.h
	${MPICXX} -o${OBJ}h5w_serial.o -c ${HDF5_PATH}/h5w_serial.cpp ${FLAGS} ${HDF5_INCL}

amuse_interface: ${OBJ}amuse_interface.o 
amuse_interface_gl: ${OBJ}amuse_interface_gl.o 

${OBJ}amuse_interface.o: ${SRC}amuse_interface.cpp ${SRC}SimpleX.cpp ${SRC}SimpleX.h ${SRC}Structs.h ${SRC}Common.h ${PLUGINS}/octree/tree_structures.cpp
	${MPICXX} -o${OBJ}amuse_interface.o ${MYINCLUDE} ${MACRODEFS} -c ${SRC}amuse_interface.cpp ${FLAGS} ${HDF5_INCL}

${OBJ}amuse_interface_gl.o: ${SRC}amuse_interface_gl.cpp ${OBJ}amuse_interface.o ${SRC}SimpleX.cpp ${SRC}SimpleX.h ${SRC}Structs.h ${SRC}Common.h ${PLUGINS}/octree/tree_structures.cpp
	${MPICXX} -o${OBJ}amuse_interface_gl.o ${MYINCLUDE} ${MACRODEFS} -c ${SRC}amuse_interface_gl.cpp ${FLAGS} ${HDF5_INCL}

clean:
	rm -rf ${OBJ}*.o ${BIN}SimpleX
