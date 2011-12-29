
# Makefile to build the worker code.  First draft was auto-generated;
# subsequent versions are hand-coded, so save a copy.  Note that
# all compiles must use MPICXX, not CXX.

# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

ALL	 = hacs64_worker hacs64_worker_cpu
OPT	  = -I$(AMUSE_DIR)/lib/stopcond -O3

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g $(OPT)
CFLAGS   += -I./src -I./src/ON_neib -fopenmp
CXXFLAGS += $(CFLAGS)  
LDFLAGS  += -L$(AMUSE_DIR)/lib/stopcond -lstopcond -lm $(MUSE_LD_FLAGS) -fopenmp

CUDA_TK  ?= /usr/local/cuda

CUDA_LIBDIRS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64
CUDA_LIBS ?= $(CUDA_LIBDIRS) -lcudart

LDFLAGS_GPU =	$(LDFLAGS) \
		$(CUDA_LIBS) \

OBJS = interface.o

CODELIB = src/libhacs64gpu.a
CODELIBCPU = src/libhacs64cpu.a

AMUSE_DIR ?= ../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: $(ALL)

save:
	@ for f in Makefile interface.py interface.cc src/Makefile; do \
		cp -p $$f $$f.save; done

$(CODELIB):
	make -C src all CXX=$(MPICXX)
	make -C src cpu CXX=$(MPICXX)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py hacs64Interface -o $@
	
interface.h: interface.py
	$(CODE_GENERATOR) --type=H  -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py hacs64Interface -o $@

hacs64_worker: worker_code.cc interface.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) worker_code.cc $(OBJS) $(CODELIB) -o $@ $(LDFLAGS_GPU) -DGPU

hacs64_worker_cpu: worker_code.cc interface.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) worker_code.cc $(OBJS) $(CODELIBCPU) -o $@  $(LDFLAGS)

%.o: %.cc
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 

clean:
	$(RM) -r -f *.so *.dSYM *.o *.mpio *.pyc worker_code.cc interface.h 
	$(RM) -f *~ hacs64_worker hacs64_worker_cpu worker_code.cc
	$(RM) -f worker_code-sockets.cc
	make -C src clean
