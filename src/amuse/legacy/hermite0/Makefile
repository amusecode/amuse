CC = mpic
CXX = mpicxx

CFLAGS ?= -g -Wall -DTOOLBOX  $(MUSE_INCLUDE_DIR) 
CXXFLAGS ?= $(CFLAGS) 

LDFLAGS = -lm $(MUSE_LD_FLAGS)

A_OBJS = hermite_code.o

AMUSE_DIR ?= ../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

SCLIBS ?= -L../../../../lib/stopcond -lstopcond
SCINC ?= -I../../../../lib/stopcond

all: worker_code

cleanall: clean 
	$(RM) muse_worker *~ worker_code

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py HermiteInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py HermiteInterface -o $@

worker_code:	worker_code.cc worker_code.h $(A_OBJS)
	mpicxx $(CXXFLAGS) $(SCINC) $@.cc $(A_OBJS) -o $@ $(SCLIBS) 

.cc.o: $<
	$(CXX) $(CXXFLAGS) $(SCINC) -c -o $@ $< 

.C.o:
	$(CC) $(CFLAGS) $(SCINC) -c -o $@ $?

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h
