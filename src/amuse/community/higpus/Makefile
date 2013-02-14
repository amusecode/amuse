

# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

# Set to empty string on purpose
# otherwhise will cause error in main makefile
# for most users (most users have no gpu).
# to easily create code in this directory run:
# make code
ALL	 = 

# all the ?= lines work as follows
# if the variable is already defined, it will keep it's pre-defined value.
# if the variable is not defined it will get the value after the ?= symbol.
MPICXX   ?= mpicxx

CUDA_INSTALL_PATH := /usr/local/cuda
MPI_INSTALL_PATH := /usr/lib64/openmpi/1.4-gcc

NVCC ?= /usr/local/cuda/bin/nvcc

#possible options : -DCHECK_ERRORS -DCHECK_TIMES -DPLUMMER -DUNROLL -DGPUCORR -DGALAXY
MYOPTS := -DUNROLL -DCHECK_ERRORS 

CUDAINCLUDE := -I$(CUDA_INSTALL_PATH)/include
CUDALIB := -L$(CUDA_INSTALL_PATH)/lib64 -lcuda -lcudart

MPIINCLUDE := -I$(MPI_INSTALL_PATH)/include
MPILIB := -L$(MPI_INSTALL_PATH)/lib

MPI_CXXFLAGS ?= $(MPIINCLUDE)
MPI_CXXLIBS ?= $(MPILIB)  -lmpi -lmpi_cxx

OPENMP_CFLAGS ?= -fopenmp

MYINCLUDE := -I./lib/ 

## for amuse, do not set a fixed LIBS line like this
## LIBS := -lcuda -lcudart -lmpi -lmpi_cxx -lm -fopenmp
## But make it depend on flags set in the config.mk file
## (note all these variables have defaults in this makefile too, for unconfigured makes)
LIBS := -lm $(OPENMP_CFLAGS) $(MPI_CXXLIBS) $(CUDALIB)

GXXOPTS := -O3 -fopenmp -Wall -W
NVCCFLAGS := -O3 --compiler-bindir=/usr/bin/$(CXX) -Xcompiler "$(GXXOPTS) $(MYOPTS)" -arch sm_13 $(CUDAINCLUDE) $(MPI_CXXFLAGS) $(MYINCLUDE)
CXXFLAGS += $(GXXOPTS) $(MYOPTS) $(CUDAINCLUDE) $(MPI_CXXFLAGS) $(MYINCLUDE)

LDFLAGS = $(LIBS) $(CUDALIB) 

OBJS = interface.o

CODELIB = src/libhigpus.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: $(ALL)

code: higpus_worker_gpu

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
	$(RM) *~ higpus_worker higpus_worker_gpu worker_code.cc
	make -C src clean
