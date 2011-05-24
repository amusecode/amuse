MPICXX ?= mpicxx
MPICC ?= mpicc

CC = $(MPICC)
CXX = $(MPICXX)

CFLAGS ?= -O3 -Wall -DTOOLBOX  $(MUSE_INCLUDE_DIR) 
CXXFLAGS ?= $(CFLAGS) 

LDFLAGS = -lm $(MUSE_LD_FLAGS)

A_OBJS = hermite_code.o

AMUSE_DIR ?= ../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

all: hermite_worker

cleanall: clean 
	$(RM) hermite_worker *~ worker_code hermite_worker_sockets

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py HermiteInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py HermiteInterface -o $@

hermite_worker:	worker_code.cc worker_code.h $(A_OBJS)
	$(MPICXX) $(CXXFLAGS) $(SCINC) $< $(A_OBJS) -o $@ $(SCLIBS) 

worker_code-sockets.cc: interface.py
	$(CODE_GENERATOR) --type=c --mode=sockets interface.py HermiteInterface -o $@

hermite0_worker_sockets: worker_code-sockets.cc worker_code.h $(CODELIB) $(A_OBJS)
	$(MPICXX) $(CXXFLAGS) $(SCINC) $(A_OBJS) $(CODELIB) -o $@ $(SCLIBS) $<

.cc.o: $<
	$(CXX) $(CXXFLAGS) $(SCINC) -c -o $@ $< 

.C.o:
	$(CC) $(CFLAGS) $(SCINC) -c -o $@ $?

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h
	rm -f hermite_worker
