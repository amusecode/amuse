# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(strip $(AMUSE_DIR)),)
    AMUSE_DIR = ../../../..
endif
-include ${AMUSE_DIR}/config.mk

ifeq ($(strip $(MPICXX)),)
    MPICXX = mpicxx
endif
ifeq ($(strip $(MPICC)),)
    MPICC = mpicc
endif

CC = $(MPICC)
CXX = $(MPICXX)

ifeq ($(strip $(CFLAGS)),)
    CFLAGS = -O3 -Wall -DTOOLBOX
endif
ifeq ($(strip $(CXXFLAGS)),)
    CXXFLAGS = $(CFLAGS) 
endif

LDFLAGS += -lm $(MUSE_LD_FLAGS)

A_OBJS = interface.o
AM_LIBS = -L$(AMUSE_DIR)/lib/amuse_mpi -lamuse_mpi
AM_CFLAGS = -I$(AMUSE_DIR)/lib/amuse_mpi
CODE_GENERATOR = $(PYTHON) $(AMUSE_DIR)/build.py

#ifeq ($(strip $(PROFLIBS)),)
#    PROFLIBS = -L$(AMUSE_DIR)/lib/ibis/src/profiling -libisprofiling -lmpi
#endif

all: symple_worker

cleanall: clean 
	$(RM) symple_worker *~ worker_code

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py sympleInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py sympleInterface -o $@

symple_worker:	worker_code.cc worker_code.h $(A_OBJS)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(AM_CFLAGS) $(LDFLAGS) $< $(A_OBJS) -o $@ $(SC_MPI_CLIBS)  $(LIBS)  $(AM_LIBS)

.cc.o: $<
	$(CXX) $(CXXFLAGS) $(SC_FLAGS) $(AM_CFLAGS) -c -o $@ $< 

.C.o:
	$(CC) $(CFLAGS) $(SC_FLAGS) $(AM_CFLAGS) -c -o $@ $?

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h
	rm -f symple_worker 
	rm -rf *.dSYM
