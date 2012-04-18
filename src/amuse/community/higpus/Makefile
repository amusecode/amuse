

# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

ALL	 = 

MPICXX   ?= mpicx

CUDA_INSTALL_PATH := /usr/local/cuda
MPI_INSTALL_PATH := /usr/lib64/openmpi/1.4-gcc

NVCC := /usr/local/cuda/bin/nvcc

#possible options : -DCHECK_ERRORS -DCHECK_TIMES -DPLUMMER -DUNROLL
MYOPTS := -DUNROLL -DCHECK_TIMES -DCHECK_ERRORS

CUDAINCLUDE := -I$(CUDA_INSTALL_PATH)/include
CUDALIB := -L$(CUDA_INSTALL_PATH)/lib64

MPIINCLUDE := -I$(MPI_INSTALL_PATH)/include
MPILIB := -L$(MPI_INSTALL_PATH)/lib

MYINCLUDE := -I./lib/ -I./src/lib

LIBS := -lcuda -lcudart -lmpi -lmpi_cxx -lm -fopenmp

GXXOPTS := -O3 -fopenmp -mno-fused-madd -m64 -Wall -W
NVCCFLAGS := -O3 --compiler-bindir=/usr/bin/$(CXX) -Xcompiler "$(GXXOPTS) $(MYOPTS)" -arch sm_13 $(CUDAINCLUDE) $(MPIINCLUDE) $(MYINCLUDE)
CXXFLAGS := $(GXXOPTS) $(MYOPTS) $(CUDAINCLUDE) $(MPIINCLUDE) $(MYINCLUDE)

LDFLAGS = $(CUDALIB) $(MPILIB) $(LIBS) 

OBJS = interface.o

CODELIB = src/libhigpus.a

AMUSE_DIR ?= ../

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: $(ALL)


$(CODELIB):
	make -C src all 

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py HiGPUsInterface -o $@
	
worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py HiGPUsInterface -o $@

higpus_worker_gpu: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) worker_code.cc $(OBJS) $(CODELIB) -o $@ $(LDFLAGS)

%.o: %.cc
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) *~ higpus_worker worker_code.cc
	make -C src clean
