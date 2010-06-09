.SUFFIXES: .cu 

CUDAPATH    := /data2/vanelteren/cuda
CUDASDKPATH := /data2/vanelteren/cuda-sdk/C
CUDAINCLUDE := -I$(CUDAPATH)/include -I$(CUDASDKPATH)/common/inc

CXX  := g++
CC   := gcc
NVCC := $(CUDAPATH)/bin/nvcc


NVCCFLAGS :=  -D_CONSOLE -O0 -g -D_DEBUG  -maxrregcount=32 $(CUDAINCLUDE)
CXXFLAGS  := -g -DNGB $(CUDAINCLUDE)
LDFLAGS   := 

OBJS = \
	sapporo.o \
	send_fetch_data.o \
	sapporoG6lib.o 

CUOBJS = host_evaluate_gravity.cu_o
TARGET = libsapporo.a

all:   $(TARGET) $(PROG)

$(TARGET): $(OBJS) $(CUOBJS)
	ar qv $@ $^
	ranlib $@

.cpp.o: 
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.cu_o:  %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	/bin/rm -rf *.o *.cu_o $(TARGET)
	/bin/rm -rf test_gravity_block test_gravity_N2ngb

