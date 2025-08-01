.SUFFIXES: .cu 

CXX  := g++
CC   := gcc
NVCC := nvcc
AR   := ar ruv
RANLIB := ranlib

CUDAPATH       := /usr/local/cuda
CUDAINCLUDE    := -I$(CUDAPATH)/include
CUDASDKINCLUDE := -I$(HOME)/NVIDIA_GPU_Computing_SDK/C/common/inc

NVCCFLAGS := -arch=sm_20 -O4 -g  $(CUDAINCLUDE) $(CUDASDKINCLUDE) -I./ -Xptxas -v,-abi=no 
# NVCCFLAGS += -maxrregcount=32
NVCCFLAGS += -Xcompiler="-Wall"
#NVCCFLAGS += -Xopencc="-O0"
#NVCCFLAGS += -G
#NVCCFLAGS := -O0 -g -D_DEBUG -deviceemu -maxrregcount=32 $(CUDAINCLUDE)


CXXFLAGS  := -O3 -g -Wall  -I./ON_neib -I./ -fopenmp
LDFLAGS   :=  -L/usr/local/cuda/lib64 -lcudart -fopenmp

OBJS = hacs64_amuse.o \
       irrf6.o \
			 regf4_cu/regf4.cu_o 

TARGET = integrate-gpu
CODELIB = libhacs64.a
all: $(CODELIB) $(TARGET) 

$(TARGET): $(CODELIB) integrate.o
	$(CXX) -o $(TARGET) integrate.o $(LDFLAGS) ./$(CODELIB) 

$(CODELIB): $(OBJS)
	/bin/rm -f $@
	$(AR) $@ $(OBJS)

.cpp.o: 
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.cu_o:  %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	/bin/rm -rf *.o regf4_cu/*.cu_o $(TARGET) $(CODELIB)

$(CODELIB): hacs64.h irrf6.h regf4.h hacs4_force.h hacs6_force.h hacs64_particle.h Scheduler.h
integrate.o: $(CODELIB)
irrf6.o: irrf6.h
regf4_cu/regf4.cu_o: regf4.h regf4_cu/cuda_pointer.h regf4_cu/dev_regf4.cu
regf4_cu/regf4.ptx: regf4.h regf4_cu/cuda_pointer.h regf4_cu/dev_regf4.cu

