# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

###############################################################
# Download and install the following libraries:
# gmp: https://gmplib.org/
# mpfr: http://www.mpfr.org/
# mpfrc++: http://www.holoborodko.com/pavel/mpfr/
# Change these paths in this Makefile and in src/makefile
GMP_LIBS   ?= -L/data1/boekholt/Packages/mpreal_local/lib/ -lgmp  # libgmp.a
MPFR_LIBS  ?= -L/data1/boekholt/Packages/mpreal_local/lib/ -lmpfr # libmpfr.a
MPFR_FLAGS ?= -I/data1/boekholt/Packages/mpreal_local/include/    # mpreal.h
###############################################################

LIBS += $(MPFR_LIBS) $(GMP_LIBS)

MPICXX ?= mpicxx
MPICC ?= mpicc

#CC = $(MPICC)
#CXX = $(MPICXX)

CXXFLAGS ?= -Wall -g -O2 

CXXFLAGS += $(MPFR_FLAGS) -std=c++0x

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CODELIB = src/libbrutus.a
OBJ = interface.o

all: brutus_worker 

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h 
	rm -f *~ brutus_worker
	make -C src clean

$(CODELIB):
	make -C src all CXXFLAGS="$(CXXFLAGS)"

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py BrutusInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py BrutusInterface -o $@

brutus_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJ)
	$(MPICXX) $(CXXFLAGS) $(MPFR_FLAGS) $(SC_FLAGS) -I./src $< $(CODELIB) $(OBJ) -o $@  -L./src -lbrutus $(SC_CLIBS) $(LIBS)

interface.o: interface.cc
	$(MPICXX) $(CXXFLAGS) $(MPFR_FLAGS) $(SC_FLAGS) -I./src -c -o $@ $< 


