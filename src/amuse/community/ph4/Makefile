
# Makefile to build the worker code.  First draft was auto-generated;
# subsequent versions are hand-coded, so save a copy.  Note that
# all compiles must use MPICXX, not CXX.

ALL	 = ph4_worker ph4_worker_gpu
OPT	  = -O4

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g $(OPT)
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CUDA_TK  ?= /usr/local/cuda
SAPPOROLIBS ?= -L/usr/local/sapporo_v1.5 -lsapporo
BOOSTLIBS ?= -L/opt/local/lib -lboost_thread-mt

LDFLAGS_GPU = $(LDFLAGS) \
		$(SAPPOROLIBS) \
		-L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64 -lcudart \
		$(BOOSTLIBS)

OBJS = interface.o

CODELIB = src/libph4.a

CODELIB_GPU = src/libph4gpu.a

CUDA_LIBS ?= -lcudart

SAPPORO_LIBDIRS ?= -L$(AMUSE_DIR)/lib/sapporo_light

AMUSE_DIR ?= ../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: $(ALL)

save:
	@ for f in Makefile interface.py interface.cc src/Makefile; do \
		cp -p $$f $$f.save; done

$(CODELIB):
	make -C src all CXX=$(MPICXX)

$(CODELIB_GPU):
	make -C src gpulib CXX=$(MPICXX)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py ph4Interface -o $@

interface.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py ph4Interface -o $@

ph4_worker: worker_code.cc interface.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) worker_code.cc $(OBJS) $(CODELIB) -o $@

ph4_worker_gpu: worker_code.cc interface.h $(CODELIB_GPU) $(OBJS)
	$(MPICXX) $(CXXFLAGS) worker_code.cc $(OBJS) $(CODELIB_GPU) -o $@ $(LDFLAGS_GPU) -DGPU

%.o: %.cc
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 

clean:
	$(RM) -r -f *.so *.dSYM *.o *.pyc worker_code.cc interface.h 
	$(RM) -f *~ ph4_worker ph4_worker_gpu worker_code.cc
	make -C src clean
