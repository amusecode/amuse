# standard amuse configuration include
# config.mk will be made after ./configure has run
# # don't run make directly in the PeTar/amuse-interface directory, this is used inside amuse.
AMUSE_DIR?=$(shell amusifier --get-amuse-dir)
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

INCLUDE  += -I src/PeTar/src -I src/SDAR/src -I src/FDPS/src 

MPIFLAGS +=  -DPARTICLE_SIMULATOR_MPI_PARALLEL -DMPICH_IGNORE_CXX_SEEKC

CXXFLAGS = ${INCLUDE} -fPIC -O2 -Wall -std=c++11  -fopenmp

# arch
CXXFLAGS += -march=core-avx2
CXXFLAGS += -D INTRINSIC_X86
CXXFLAGS += -D USE_SIMD
CXXFLAGS += -D DIV_FIX

CXXFLAGS += -D PARTICLE_SIMULATOR_THREAD_PARALLEL 
# MPI
CXXFLAGS += -D PARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -D MPICH_IGNORE_CXX_SEEKC

CXXFLAGS += -D SOFT_PERT -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D CLUSTER_VELOCITY
CXXFLAGS += -D USE_QUAD
CXXFLAGS += -D STELLAR_EVOLUTION

CXXFLAGS += -D PROFILE
CXXFLAGS += -D HARD_CHECK_ENERGY
CXXFLAGS += -D TIDAL_TENSOR_3RD

# debug
#CXXFLAGS +=  -D PETAR_DEBUG -D AR_DEBUG -D AR_DEBUG_DUMP -D AR_WARN -D HARD_DEBUG -D HARD_DUMP -D CLUSTER_DEBUG -D ARTIFICIAL_PARTICLE_DEBUG 
CXXFLAGS += -D INTERFACE_DEBUG_PRINT
CXXFLAGS += -D INTERFACE_DEBUG

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o 

FILELIST=interface.cc interface.py interface.h

CODELIB = 

all: petar_worker __init__.py

clean:
	$(RM) -rf __pycache__
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ petar_worker worker_code.cc
#	make -C src clean

distclean: clean
	make -C src $(FILELIST)

#$(CODELIB):
#	make -C src all

DOWNLOAD_FROM_WEB = $(PYTHON) ./download.py

download:
	$(RM) -Rf .pc
	$(RM) -Rf src
	mkdir src
	$(DOWNLOAD_FROM_WEB)
	ln -sf src/PeTar/amuse-interface/interface.py ./
	ln -sf src/PeTar/amuse-interface/interface.h ./
	ln -sf src/PeTar/amuse-interface/interface.cc ./



__init__.py:
	touch $@

$(FILELIST): |download

test_interface: $(OBJS) test_interface.o
	$(MPICXX) $(CXXFLAGS) $(OBJS) test_interface.o -o $@

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py petarInterface -o $@

worker_code.h: interface.py 
	$(CODE_GENERATOR) --type=H -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py petarInterface -o $@

petar_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS) interface.h
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(AM_CFLAGS) $(LDFLAGS) $< $(OBJS) $(CODELIB) -o $@ $(SC_MPI_CLIBS) $(LIBS) $(AM_LIBS)

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) -c -o $@ $< 
