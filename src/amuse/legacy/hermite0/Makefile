CC = mpic
CXX = mpicxx

CFLAGS ?= -O3 -Wall -DTOOLBOX  $(MUSE_INCLUDE_DIR) 
CXXFLAGS ?= $(CFLAGS) 

LDFLAGS = -lm $(MUSE_LD_FLAGS)

A_OBJS = hermite_code.o

AMUSE_DIR ?= ../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

all: worker_code

cleanall: clean 
	$(RM) muse_worker *~ worker_code

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py HermiteInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H -i amuse.support.legacy.stopping_conditions.StoppingConditionInterface interface.py HermiteInterface -o $@

worker_code:	worker_code.cc worker_code.h $(A_OBJS)
	mpicxx $(CXXFLAGS) bla $(SCINC) $@.cc $(A_OBJS) -o $@ $(SCLIBS) 

.cc.o: $<
	$(CXX) $(CXXFLAGS) $(SCINC) -c -o $@ $< 

.C.o:
	$(CC) $(CFLAGS) $(SCINC) -c -o $@ $?

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h
