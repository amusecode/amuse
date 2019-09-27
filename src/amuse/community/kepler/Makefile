# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

# Makefile to build the worker code.  First draft was auto-generated;
# subsequent versions are hand-coded, so save a copy.  Note that
# all compiles must use MPICXX, not CXX.

ALL	 = kepler_worker
OPT	  = # -O4 don't use O4 it is not defined on normal gcc, but has special meaning on LLVM/clang gcc

AM_LIBS = -L$(AMUSE_DIR)/lib/amuse_mpi -lamuse_mpi
AM_CFLAGS = -I$(AMUSE_DIR)/lib/amuse_mpi

MPICXX   ?= mpicxx

CFLAGS   += $(OPT) -Isrc
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libkepler.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: $(ALL)

$(CODELIB):
	make -C src all CXX='$(MPICXX)'

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py KeplerInterface -o $@

interface.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py KeplerInterface -o $@

kepler_worker: worker_code.cc interface.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) worker_code.cc $(OBJS) $(CODELIB) -o $@ $(LIBS)

kepler_worker_cython: kepler_cython.so
	$(CODE_GENERATOR) --type=cython -m script -x amuse.community.kepler.interface KeplerInterface -o $@ --cython-import kepler_cython
	
kepler_cython.so: kepler_cython.o $(CODELIB) $(OBJS)
	$(MPICXX) -shared $(CXXFLAGS) $(PYTHONDEV_LDFLAGS) $(AM_CFLAGS) $(SC_FLAGS) $(LDFLAGS) kepler_cython.o  -o $@ $(SC_CLIBS) $(AM_LIBS) $(OBJS) $(CODELIB)

kepler_cython.o: kepler_cython.cc
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(AM_CFLAGS) $(PYTHONDEV_CFLAGS) -c -o $@ $< 

kepler_cython.pyx: interface.py interface.h
	$(CODE_GENERATOR) --type=cython -m mpi interface.py KeplerInterface -o $@

kepler_cython.cc: kepler_cython.pyx
	$(CYTHON) $< -o $@



%.o: %.cc
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 

clean:
	$(RM) -r -f *.so *.dSYM *.o *.pyc worker_code.cc interface.h 
	$(RM) -f *~ kepler_worker worker_code.cc kepler_worker_sockets worker_code-sockets.cc
	$(RM) -f kepler_worker_cython*
	make -C src clean

distclean: clean
