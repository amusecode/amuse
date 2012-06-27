# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

# Makefile to build the worker code.  First draft was auto-generated;
# subsequent versions are hand-coded, so save a copy.  Note that
# all compiles must use MPICXX, not CXX.

ALL	 = smallN_worker
OPT	  = -I$(AMUSE_DIR)/lib/stopcond

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g $(OPT) -I../../ph4/src -I../ph4/src
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libsmallN.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: $(ALL)

$(CODELIB):
	make -C src all CXX=$(MPICXX)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py SmallNInterface -o $@

interface.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py SmallNInterface -o $@

smallN_worker: worker_code.cc interface.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) worker_code.cc $(OBJS) $(CODELIB) -o $@ $(LDFLAGS)

%.o: %.cc
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 

clean:
	$(RM) -r -f *.so *.dSYM *.o *.pyc worker_code.cc interface.h 
	$(RM) -f *~ smallN_worker worker_code.cc
	make -C src clean
