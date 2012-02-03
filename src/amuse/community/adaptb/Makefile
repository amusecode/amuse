# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

LDIR = /home/boekholt/Packages/
LLDIR = /home/boekholt/Packages/mpfr_c++/
LIBS = -L$(LDIR)  -lmpfr -lgmp

MPICXX   ?= mpicxx
CXXFLAGS ?= -Wall -g -O2
#CXXFLAGS += -I${LLDIR}

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcond
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CODELIB = src/libadaptb.a

OBJS = interface.o

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

adaptb_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) -I${LLDIR} $(SCINC) -I./src $^ -o $@  -L./src -ladaptb $(SCLIBS) $(LIBS)

#adaptb_worker: worker_code.o $(CODELIB) $(OBJS)
#	$(MPICXX) $(CXXFLAGS) $(SCINC) $< $(CODELIB) $(OBJS) -o $@ $(SCLIBS) $(LIBS) -L./src -ladaptb

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) -I${LLDIR} $(SCINC) -I./src -c -o $@ $< 

#worker_code.o: worker_code.cc worker_code.h 
#	$(MPICXX) $(CXXFLAGS) -I${LLDIR} $(SCINC) -I./src -c -o $@ $< 

