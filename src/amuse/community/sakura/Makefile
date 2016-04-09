# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx
MPICC ?= mpicc

CC = $(MPICC)
CXX = $(MPICXX)

CXXFLAGS ?= -Wall -g -O2

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CODELIB = src/libsakura.a
OBJ = interface.o

all: sakura_worker 

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h 
	rm -f *~ sakura_worker
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py SakuraInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py SakuraInterface -o $@

sakura_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJ)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(LDFLAGS) -I./src $< $(CODELIB) $(OBJ) -o $@  -L./src -lsakura $(SC_CLIBS)

interface.o: interface.cc
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) -I./src -c -o $@ $< 

