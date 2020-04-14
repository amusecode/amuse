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

CFLAGS ?= -O3 -Wall -DTOOLBOX
CXXFLAGS ?= $(CFLAGS) 

LDFLAGS += -lm $(MUSE_LD_FLAGS)

A_OBJS = interface.o
CODE_GENERATOR ?= $(PYTHON) $(AMUSE_DIR)/build.py

#PROFLIBS ?= -L$(AMUSE_DIR)/lib/ibis/src/profiling -libisprofiling -lmpi

all: symple_worker

cleanall: clean 
	$(RM) symple_worker *~ worker_code

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py sympleInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py sympleInterface -o $@

symple_worker:	worker_code.cc worker_code.h $(A_OBJS)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(LDFLAGS) $< $(A_OBJS) -o $@ $(SC_MPI_CLIBS)  $(LIBS)

.cc.o: $<
	$(CXX) $(CXXFLAGS) $(SC_FLAGS) $(AM_CFLAGS) -c -o $@ $< 

.C.o:
	$(CC) $(CFLAGS) $(SC_FLAGS) $(AM_CFLAGS) -c -o $@ $?

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h
	rm -f symple_worker 
	rm -rf *.dSYM
