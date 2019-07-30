# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

# Makefile to build the worker code.  First draft was auto-generated;
# subsequent versions are hand-coded, so save a copy.  Note that
# all compiles must use MPICXX, not CXX.

ALL	 = keplerorbiters_worker
OPT	  = $(OPENMP_CFLAGS) # -O4 don't use O4 it is not defined on normal gcc, but has special meaning on LLVM/clang gcc
#OPT	  =

MPICXX   ?= mpicxx

CFLAGS   += $(OPT)
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

keplerorbiters_worker: worker_code.cc interface.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) worker_code.cc $(OBJS) $(CODELIB) -o $@ $(LIBS)

%.o: %.cc
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 

clean:
	$(RM) -r -f *.so *.dSYM *.o *.pyc worker_code.cc interface.h 
	$(RM) -f *~ kepler_worker keplerorbiters_worker worker_code.cc worker_code-sockets.cc
	make -C src clean
