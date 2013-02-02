# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../..
-include $(AMUSE_DIR)/config.mk

.SUFFIXES: .cpp .cu .o

CODELIB = libhigpus.a

MPICXX    ?= mpicxx
MPICC     ?= mpicc
CXX	  ?= mpicxx

CUDA_INSTALL_PATH := /usr/local/cuda
MPI_INSTALL_PATH := /usr/lib64/openmpi/1.4-gcc

NVCC := /usr/local/cuda/bin/nvcc

#possible options : -DCHECK_ERRORS -DCHECK_TIMES -DPLUMMER -DUNROLL
MYOPTS := -DUNROLL -DCHECK_ERRORS

CUDAINCLUDE := -I$(CUDA_INSTALL_PATH)/include
CUDALIB := -L$(CUDA_INSTALL_PATH)/lib64

MPIINCLUDE := -I$(MPI_INSTALL_PATH)/include
MPILIB := -L$(MPI_INSTALL_PATH)/lib

MPI_CXXFLAGS = $(MPIINCLUDE)
MPI_CXXLIBS = $(MPILIB) -lmpi -lmpi_cxx

MYINCLUDE := -I./lib/

LIBS := -lcuda -lcudart -lmpi -lmpi_cxx -lm -fopenmp

GXXOPTS := -O3 -fopenmp -Wall -W
NVCCFLAGS := -O3 --compiler-bindir=/usr/bin/$(CXX) -Xcompiler "$(GXXOPTS) $(MYOPTS)" -arch sm_13 $(CUDAINCLUDE) $(MPI_CXXFLAGS) $(MYINCLUDE)
CXXFLAGS := $(GXXOPTS) $(MYOPTS) $(CUDAINCLUDE) $(MPI_CXXFLAGS) $(MYINCLUDE)

CUDAOBJS = 	./src/cuda/main.cu.o \
                ./src/cuda/cpu_func.cu.o \
                ./src/cuda/corrector.cu.o \
                ./src/cuda/cudainit.cu.o \
                ./src/cuda/hermite6.cu.o \
                ./src/cuda/warmstart.cu.o \
                ./src/cuda/selectactive.cu.o \
                ./src/cuda/initblocks.cu.o \
                ./src/cuda/allreduce.cu.o \
                ./src/cxx/utilis.cu.o \
                ./src/cuda/energy.cu.o \
                ./src/cuda/cuda_kernels.cu.o

CXXOBJS = ./src/cxx/types.cpp.o

CODEOBJS = $(CUDAOBJS) $(CXXOBJS)

TARGET = ./exec/HiGPUs.x

AR = ar ruv
RANLIB = ranlib
RM = rm

LINKLINE = $(CXX) -o $(TARGET) $(CODEOBJS) $(CUDALIB) $(LIBS) $(MPI_CXXLIBS)

all: $(CODELIB)

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@


%.cu.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

exec: $(TARGET)

$(TARGET): $(CUDAOBJS) $(CXXOBJS) Makefile 
	$(LINKLINE)

clean:
	/bin/rm -rf ./src/cuda/*.o ./src/cxx/*.o ./exec/HiGPUs.x ./exec/nohup.out ./exec/times.dat
	/bin/rm -rf ./exec/data/energy.dat ./exec/data/Blocks.dat ./exec/data/HiGPUslog.dat ./exec/data/*memory*
	/bin/rm -rf ./*.a
