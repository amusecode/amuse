.SUFFIXES: .cu 

CXX  := g++
CC   := gcc
NVCC ?= nvcc
AR   := ar ruv
RANLIB := ranlib

CUDAPATH       := /usr/local/cuda
CUDAINCLUDE    := -I$(CUDAPATH)/include
#CUDASDKINCLUDE := -I$(CUDASDK)/common/inc

#NVCCFLAGS := -arch=sm_20 -O4 -g  $(CUDAINCLUDE) $(CUDASDKINCLUDE) -I./ -Xptxas -v,-abi=no 
NVCCFLAGS := -arch=sm_20 -O4 -g  $(CUDAINCLUDE)  -I./ -Xptxas -v,-abi=no 
# NVCCFLAGS += -maxrregcount=32
NVCCFLAGS += -Xcompiler="-Wall"
#NVCCFLAGS += -Xopencc="-O0"
#NVCCFLAGS += -G
#NVCCFLAGS := -O0 -g -D_DEBUG -deviceemu -maxrregcount=32 $(CUDAINCLUDE)
CUDA_LIBS = -lcudart  -L $(CUDAPATH)/lib -L $(CUDAPATH)/lib64

OPENMP_CFLAGS ?= -fopenmp

CXXFLAGS  +=  -I./ON_neib -I./  $(CUDAINCLUDE)
LDFLAGS   +=  $(OPENMP_CFLAGS)
LDGPUGLAGS := $(LDFLAGS) $(CUDA_LIBS)

OBJS = hacs64_amuse.o \
       irrf6.o \
			 regf4_cu/regf4.cu_o 

OBJSCPU = hacs64_amuse.o \
       irrf6.o \
			 regf4.o 

TARGET = integrate-gpu
CODELIB = libhacs64gpu.a

TARGETCPU  = integrate-cpu
CODELIBCPU = libhacs64cpu.a

all: $(CODELIB) $(TARGET) 

cpu: $(CODELIBCPU) $(TARGETCPU)

$(TARGET): $(CODELIB) integrate.o
	$(CXX) -o $(TARGET) integrate.o $(LDFLAGS)  $(CUDA_LIBS) ./$(CODELIB) 

$(CODELIB): $(OBJS)
	/bin/rm -f $@
	$(AR) $@ $(OBJS)

$(TARGETCPU): $(CODELIBCPU) integrate.o
	$(CXX) -o $(TARGETCPU) integrate.o $(LDFLAGS) ./$(CODELIBCPU) 

$(CODELIBCPU): $(OBJSCPU)
	/bin/rm -f $@
	$(AR) $@ $(OBJSCPU)

.cpp.o: 
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.cu_o:  %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	/bin/rm -rf *.o regf4_cu/*.cu_o $(TARGET) $(CODELIB) $(TARGETCPU) $(CODELIBCPU)

$(CODELIB): hacs64.h irrf6.h regf4.h hacs4_force.h hacs6_force.h hacs64_particle.h Scheduler.h kdtree.h
$(OBJS): hacs64.h irrf6.h regf4.h hacs4_force.h hacs6_force.h hacs64_particle.h Scheduler.h kdtree.h
integrate.o:
irrf6.o: irrf6.h
regf4_cu/regf4.cu_o: regf4.h regf4_cu/cuda_pointer.h regf4_cu/dev_regf4.cu
regf4_cu/regf4.ptx: regf4.h regf4_cu/cuda_pointer.h regf4_cu/dev_regf4.cu

