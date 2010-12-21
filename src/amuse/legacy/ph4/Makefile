
# Makefile to build the worker code.  First draft was auto-generated;
# subsequent versions are hand-coded, so save a copy.  Note that
# all compiles must use MPICXX, not CXX.

MPICXX   ?= openmpicxx

CFLAGS   += -Wall -g -DGPU
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS) \
		-L/usr/local/sapporo_v1.5 -lsapporo \
		-L/usr/local/cuda/lib -lcuda -lcudart \
		-L/opt/local/lib -lboost_thread-mt

OBJS = interface.o

CODELIB = src/libph4.a

AMUSE_DIR ?= ../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: ph4_worker

save:
	@ for f in Makefile interface.py interface.cc src/Makefile; do \
		cp -p $$f $$f.save; done

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py ph4Interface -o $@

interface.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py ph4Interface -o $@

ph4_worker: worker_code.cc interface.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) worker_code.cc $(OBJS) $(CODELIB) -o $@

.cc.o:	$<
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 

clean:
	$(RM) -r -f *.so *.dSYM *.o *.pyc worker_code.cc interface.h 
	$(RM) *~ ph4_worker worker_code.cc
	make -C src clean
