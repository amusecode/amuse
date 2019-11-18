# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=$(realpath ../../../..)
-include ${AMUSE_DIR}/config.mk

# Makefile to build the worker code.  First draft was auto-generated;
# subsequent versions are hand-coded, so save a copy.  Note that
# all compiles must use MPICXX, not CXX.

ALL	 = smallN_worker
OPT	  = -I$(AMUSE_DIR)/lib/stopcond

MPICXX   ?= mpicxx
SC_FLAGS ?= -I$(AMUSE_DIR)/lib/stopcond
SC_CLIBS  ?=  -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi 

CFLAGS   += -Wall -g $(OPT) -Isrc
CXXFLAGS += $(CFLAGS) $(SC_FLAGS)
LDFLAGS  += $(SC_CLIBS) -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

AM_LIBS = -L$(AMUSE_DIR)/lib/amuse_mpi -lamuse_mpi
AM_CFLAGS = -I$(AMUSE_DIR)/lib/amuse_mpi

CODELIB = src/libsmallN.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: $(ALL)

$(CODELIB):
	make -C src all CXX='$(MPICXX)'

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py SmallNInterface -o $@

interface.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py -i amuse.support.codes.stopping_conditions.StoppingConditionInterface  SmallNInterface -o $@

smallN_worker: worker_code.cc interface.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) worker_code.cc $(OBJS) $(CODELIB) -o $@ $(LDFLAGS) $(LIBS)


smalln_worker_cython: smalln_cython.so
	$(CODE_GENERATOR) --type=cython -m script -x amuse.community.smalln.interface SmallNInterface -o $@ --cython-import smalln_cython
	
smalln_cython.so: smalln_cython.o interface.o interface.h  $(CODELIB) $(A_OBJS)
	$(MPICXX) -shared $(CFLAGS) $(PYTHONDEV_LDFLAGS) $(AM_CFLAGS) $(SC_FLAGS) $(LDFLAGS) smalln_cython.o interface.o $(A_OBJS) $(CODELIB) -o $@ $(SC_CLIBS) $(AM_LIBS) 

smalln_cython.o: smalln_cython.c interface.h
	$(MPICC) $(CFLAGS) $(SC_FLAGS) $(AM_CFLAGS) $(PYTHONDEV_CFLAGS) -c -o $@ $< 

smalln_cython.c: smalln_cython.pyx
	$(CYTHON) $< -o $@
	
smalln_cython.pyx: interface.py
	$(CODE_GENERATOR) --type=cython -m mpi interface.py SmallNInterface -o $@

	

%.o: %.cc
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 

clean:
	$(RM) -r -f *.so *.dSYM *.o *.pyc worker_code.cc interface.h 
	$(RM) -f *~ smallN_worker worker_code.cc
	$(RM) -f smalln_worker_cython smalln_cython.so smalln_cython.o smalln_cython.c smalln_cython.pyx
	make -C src clean

distclean: clean
