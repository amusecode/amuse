
# Makefile to build the worker code.  First draft was auto-generated;
# subsequent versions are hand-coded, so save a copy.  Note that
# all compiles must use MPICXX, not CXX.

AMUSE_DIR?=../../../..
ALL	 = ph4_worker ph4_worker_gpu
OPT	  = -O4 -I$(AMUSE_DIR)/lib/stopcond

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g $(OPT)
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -L$(AMUSE_DIR)/lib/stopcond -lstopcond -lm $(MUSE_LD_FLAGS)

CUDA_TK  ?= /usr/local/cuda
SAPPOROLIBS ?= -L../../../../lib/sapporo_light/ -lsapporo

# Boost is needed by sapporo.

BOOSTLIBS ?= -L/opt/local/lib -lboost_thread-mt

CUDA_LIBDIRS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64
CUDA_LIBS ?= -lcudart

LDFLAGS_GPU =	$(LDFLAGS) \
		$(SAPPOROLIBS) \
		$(CUDA_LIBDIRS) $(CUDA_LIBS) \
		$(BOOSTLIBS)

OBJS = interface.o

CODELIB = src/libph4.a

CODELIB_GPU = src/libph4gpu.a

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
	$(CODE_GENERATOR) --type=H  -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py ph4Interface -o $@

ph4_worker: worker_code.cc interface.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) worker_code.cc $(OBJS) $(CODELIB) -o $@ $(LDFLAGS)

ph4_worker_gpu: worker_code.cc interface.h $(CODELIB_GPU) $(OBJS)
	$(MPICXX) $(CXXFLAGS) worker_code.cc $(OBJS) $(CODELIB_GPU) -o $@ $(LDFLAGS_GPU) -DGPU

%.o: %.cc
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 

clean:
	$(RM) -r -f *.so *.dSYM *.o *.pyc worker_code.cc interface.h 
	$(RM) -f *~ ph4_worker ph4_worker_gpu worker_code.cc
	make -C src clean
