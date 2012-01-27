# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX   ?= mpicxx
CXXFLAGS ?= -Wall -g -O2

SAPPORO_LIBS ?= -L$(AMUSE_DIR)/lib/sapporo_light -lsapporo
CUDA_LIBS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64 -lcudart

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcond
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

CODE_GENERATOR = $(AMUSE_DIR)/build.py


CODELIB = src/libmi6.a
CODELIB_GPU = src/libmi6_gpu.a

OBJS = interface.o

all: mi6_worker 

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h 
	rm -f *~ mi6_worker mi6_worker_gpu
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
	$(MPICXX) $(CXXFLAGS) $(SCINC) -I./src  $< $(OBJS) $(CODELIB) -o $@ -L./src -lmi6 $(SCLIBS)

mi6_worker_gpu: worker_code.cc worker_code.h $(CODELIB_GPU) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(SCINC) -I./src  $< $(OBJS) $(CODELIB_GPU) -o $@ -L./src -lmi6_gpu $(SCLIBS) $(SAPPORO_LIBS) $(CUDA_LIBS)

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) $(SCINC) -I./src -c -o $@ $< 
