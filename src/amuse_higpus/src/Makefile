# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../..
-include $(AMUSE_DIR)/config.mk

.SUFFIXES: .cpp .cu .o

CODELIB = libhigpus.a

MPICXX    ?= mpicxx
MPICC     ?= mpicc
CXX	  ?= mpicxx

NVCC ?= /usr/local/cuda/bin/nvcc

#possible options : -DCHECK_ERRORS -DCHECK_TIMES -DPLUMMER -DUNROLL
MYOPTS := -DUNROLL -DCHECK_ERRORS

CUDAINCLUDE ?= -I$(CUDA_TK)/include

MYINCLUDE := -I./lib/

LIBS := -lm -fopenmp

GXXOPTS := -O3 -fopenmp -Wall -W
NVCCFLAGS := -O3 -Xcompiler "$(GXXOPTS) $(MYOPTS) $(MPI_FCFLAGS) $(MPI_CFLAGS)" $(CUDAINCLUDE) $(MYINCLUDE) $(NVCC_FLAGS)
CXXFLAGS += $(GXXOPTS) $(MYOPTS) $(CUDAINCLUDE) $(MYINCLUDE) $(MPI_CLIBS) $(MPI_CFLAGS)

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

LINKLINE = $(CXX) -o $(TARGET) $(CODEOBJS) $(CUDA_LIBS) $(LIBS) $(MPI_CLIBS)

all: $(CODELIB)

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@


%.cu.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) $(MPI_FCFLAGS) -c $< -o $@

exec: $(TARGET)

$(TARGET): $(CUDAOBJS) $(CXXOBJS) Makefile 
	$(LINKLINE)

clean:
	/bin/rm -rf ./src/cuda/*.o ./src/cxx/*.o ./exec/HiGPUs.x ./exec/nohup.out ./exec/times.dat
	/bin/rm -rf ./exec/data/energy.dat ./exec/data/Blocks.dat ./exec/data/HiGPUslog.dat ./exec/data/*memory*
	/bin/rm -rf ./*.a
