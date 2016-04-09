# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX   ?= mpicxx
CXX = $(MPICXX)
CC = $(MPICC)

CFLAGS   += -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libsei.a
CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: sei_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ sei_worker worker_code.cc
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py SeiInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py SeiInterface -o $@

sei_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< $(OBJS) $(CODELIB) -o $@ -lm $(LIBS)
	
.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
