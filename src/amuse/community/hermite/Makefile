# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx
MPICC ?= mpicc

CC = $(MPICC)
CXX = $(MPICXX)

CFLAGS ?= -O3 -Wall -DTOOLBOX  $(MUSE_INCLUDE_DIR) 
CXXFLAGS ?= $(CFLAGS) 

LDFLAGS += -lm $(MUSE_LD_FLAGS)

A_OBJS = interface.o
ifeq ($(MPI_ENABLED),yes)
AM_LIBS = -L$(AMUSE_DIR)/lib/amuse_mpi -lamuse_mpi
AM_CFLAGS = -I$(AMUSE_DIR)/lib/amuse_mpi
endif
CODE_GENERATOR ?= $(PYTHON) $(AMUSE_DIR)/build.py

#PROFLIBS ?= -L$(AMUSE_DIR)/lib/ibis/src/profiling -libisprofiling -lmpi

all: hermite_worker

cleanall: clean 
	$(RM) hermite_worker *~ worker_code hermite_worker_sockets

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py HermiteInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py HermiteInterface -o $@

hermite_worker:	worker_code.cc worker_code.h $(A_OBJS)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(AM_CFLAGS) $(LDFLAGS) $< $(A_OBJS) -o $@ $(SC_MPI_CLIBS)  $(LIBS) $(AM_LIBS)

.cc.o: $<
	$(CXX) $(CXXFLAGS) $(SC_FLAGS) $(AM_CFLAGS) -c -o $@ $< 

.C.o:
	$(CC) $(CFLAGS) $(SC_FLAGS) $(AM_CFLAGS) -c -o $@ $?

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h
	rm -f hermite_worker hermite_worker worker_code-sockets.cc hermite_worker_sockets
	rm -rf *.dSYM
