.SUFFIXES: .cu 

CUDA_SDK ?= /disks/koppoel1/CUDA23/cudasdk/C
CUDA_TK  ?= /disks/koppoel1/CUDA23/cuda

CUDAPATH    ?= $(CUDA_TK)
CUDASDKPATH ?= $(CUDA_SDK)
CUDAINCLUDE ?= -I$(CUDAPATH)/include -I$(CUDASDKPATH)/common/inc

CXX  ?= g++
CC   ?= gcc
NVCC ?= $(CUDAPATH)/bin/nvcc


NVCCFLAGS :=  -D_CONSOLE -O0 -g -D_DEBUG  -maxrregcount=32 $(CUDAINCLUDE) -m64
CXXFLAGS  := -g -DNGB $(CUDAINCLUDE) -m64
LDFLAGS   := 

OBJS = \
	sapporo.o \
	send_fetch_data.o \
	sapporoG6lib.o 

CUOBJS = host_evaluate_gravity.cu_o
TARGET = libsapporo.a

all:  $(CUDA_TK) $(CUDA_SDK) $(TARGET) $(PROG)

$(CUDA_TK):
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "------------"
	@echo "CUDA_TK variable is not set to a valid path,"
	@echo "please set the CUDA_TK variable to the directory"
	@echo "where you installed CUDA"
	@echo "the CUDA_TK directory must contain a bin directory with the <nvcc> executable"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error


$(CUDA_SDK):
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "CUDA_SDK variable is not set to a valid path,"
	@echo "please set the CUDA_SDK variable to the directory"
	@echo "where you installed the CUDA SDK"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error


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

