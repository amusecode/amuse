# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) -O3 -ffast-math -funroll-loops
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

# Simulation control macros 
# CXXFLAGS += -DENABLE_VARIABLE_SMOOTHING_LENGTH # This does not work yet with AMUSE
CXXFLAGS += -DUSE_ENTROPY
CXXFLAGS += -DUSE_BALSARA_SWITCH
#CXXFLAGS += -DUSE_PRESCR_OF_THOMAS_COUCHMAN_1992
CXXFLAGS += -DISOTHERMAL_EOS

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODEDIR = src/FDPS

INC = -I$(CODEDIR)/src
# INC += -I$(CODEDIR)/sample/c++/nbody+sph

DOWNLOAD_FROM_WEB = $(PYTHON) ./download_http.py

CODE_GENERATOR = $(AMUSE_DIR)/build.py

OBJS = src/physical_constants.o src/mathematical_constants.o interface.o

all: fdpssn_worker

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) -f src/*.o
	$(RM) *~ fdpssn_worker worker_code.cc
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py P3TARInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py P3TARInterface -o $@

fdpssn_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(INC) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@

.cc.o: $<
	$(MPICXX) $(INC) $(CXXFLAGS) -c -o $@ $< 
