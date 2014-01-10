# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

GMP_LIBS ?= -L/home/boekholt/Packages/ -lgmp
MPFR_LIBS ?= -L/home/boekholt/Packages/ -lmpfr
MPFR_FLAGS ?= -I/home/boekholt/Packages/
LIBS += $(MPFR_LIBS) $(GMP_LIBS)

MPICXX   ?= mpicxx
CXX ?= g++

CXXFLAGS ?= -Wall -g -O2

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CODELIB = src/libadaptb.a
OBJ = interface.o

all: adaptb_worker 

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h 
	rm -f *~ adaptb_worker
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py AdaptbInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py AdaptbInterface -o $@

adaptb_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJ)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) -I./src $< $(CODELIB) $(OBJ) -o $@  -L./src -ladaptb $(SC_CLIBS) $(LIBS)

interface.o: interface.cc
	$(CXX) $(CXXFLAGS) $(MPFR_FLAGS)  $(SC_FLAGS) -I./src -c -o $@ $< 

