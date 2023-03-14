# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk

GMP_LIBS ?= -lgmp
MPFR_LIBS ?= -lmpfr
LIBS += $(MPFR_LIBS) $(GMP_LIBS)

MPICXX   ?= mpicxx
CXX ?= g++

CXXFLAGS ?= -Wall -g -O2
CXXFLAGS += -std=c++11

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
	$(CODE_GENERATOR) --type=h -i amuse.community.interface.stopping_conditions.StoppingConditionInterface interface.py AdaptbInterface -o $@

adaptb_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJ)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(LDFLAGS) -I./src $< $(CODELIB) $(OBJ) -o $@  -L./src -ladaptb $(SC_CLIBS) $(LIBS)

interface.o: interface.cc
	$(CXX) $(CXXFLAGS) $(MPFR_FLAGS)  $(SC_FLAGS) -I./src -c -o $@ $< 

