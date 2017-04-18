# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk
GPUARCH ?= 35
CUDAFLAGS += -arch=sm_$(GPUARCH)
######### ugly!!!! should be in the amuse base make file


MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libetics.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: etics_worker

clean:
	$(RM) *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ etics_worker worker_code.cc
	make -C src clean

$(CODELIB):
	make -C src libetics.a

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py EticsInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py EticsInterface -o $@

etics_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) -L/usr/local/cuda/lib64 $< $(OBJS) $(CODELIB) -o $@ -lcudart
	echo the above compilation directive is not nice

.SUFFIXES: .cu .o

.cu.o: $<
	$(NVCC) -Xcompiler="$(CXXFLAGS)" -c -o $@ $<


# .cc.o: $<
# 	$(CXX) $(CXXFLAGS) -c -o $@ $<
