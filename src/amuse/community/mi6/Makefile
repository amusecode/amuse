# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX   ?= mpicxx
CXXFLAGS ?= -Wall -g -O2

SAPPORO_LIBS ?= -L$(AMUSE_DIR)/lib/sapporo_light -lsapporo
CUDA_LIBS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64 -lcudart

CODE_GENERATOR = $(AMUSE_DIR)/build.py


CODELIB = src/libmi6.a
CODELIB_GPU = src/libmi6_gpu.a

ifneq (,$(findstring g++, $(notdir $(CXX))))
CXXFLAGS += -O3 -ffast-math -funroll-loops -fomit-frame-pointer
endif


ifeq ($(findstring gcc, $(notdir $(CXX))), gcc)
CXXFLAGS += -O3 -ffast-math -funroll-loops -fomit-frame-pointer
endif

ifeq ($(findstring xlC, $(notdir $(CXX))), xlC)
CXXFLAGS += -qunroll  
endif

ifeq ($(findstring ifort, $(notdir $(CXX))), ifort)
CXXFLAGS += -O3 -ffast-math -funroll-loops -fomit-frame-pointer
endif


OBJS = interface.o

all: mi6_worker 

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h 
	rm -f *~ mi6_worker mi6_worker_gpu
	rm -f worker_code-sockets.cc mi6_worker_sockets mi6_worker_gpu_sockets
	make -C src clean

$(CODELIB):
	make -C src emulib

$(CODELIB_GPU):
	make -C src gpulib

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py MI6Interface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py MI6Interface -o $@

mi6_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) -I./src  $< $(OBJS) $(CODELIB) -o $@ -L./src -lmi6 $(SC_CLIBS)  $(LIBS)

mi6_worker_gpu: worker_code.cc worker_code.h $(CODELIB_GPU) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) -I./src  $< $(OBJS) $(CODELIB_GPU) -o $@ -L./src -lmi6_gpu $(SC_CLIBS) $(SAPPORO_LIBS) $(CUDA_LIBS)  $(LIBS)

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) -I./src -c -o $@ $< 
