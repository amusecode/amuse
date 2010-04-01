$(warning debugging this makefile)
CODEPATH = src/octgrav_v1.7d
CODELIB = $(CODEPATH)/liboctgrav.a 

CUDASDK = /data1/cuda_sdk
CUDATK  = /data1/cuda

MUSE_INCLUDE_DIR = -I$(CODEPATH) -I$(CUDASDK)/common/inc -I$(CUDATK)/include 
$(warning $(MAKE))
CXXFLAGS+= $(MUSE_INCLUDE_DIR) 


OBJS = integrator.o octgrav_code.o 

AMUSE_DIR?=../../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py

LIBS += -L$(CUDATK)/lib -lcudart

all: compile worker_code

compile: octgrav 

octgrav: 
	$(MAKE) -C $(CODEPATH)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py OctgravInterface -o $@ 

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py OctgravInterface -o $@

worker_code: worker_code.cc worker_code.h $(OBJS)
	mpicxx $(CXXFLAGS) $@.cc $(OBJS) $(LIBS) $(CODELIB) -o $@

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) worker_code *~

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

.C.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<


.PHONY: octgrav
