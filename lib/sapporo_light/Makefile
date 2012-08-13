# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../..
-include ${AMUSE_DIR}/config.mk

.SUFFIXES: .cu 

CUDA_TK  ?= /disks/koppoel1/CUDA23/cuda

CUDAPATH    ?= $(CUDA_TK)
CUDAINCLUDE ?= -I$(CUDAPATH)/include 

CXX  ?= g++
CC   ?= gcc
NVCC ?= $(CUDAPATH)/bin/nvcc


NVCCFLAGS :=  $(NVCC_FLAGS) -D_CONSOLE -O0 -g -D_DEBUG  -maxrregcount=32 $(CUDAINCLUDE)
CXXFLAGS  := -g -DNGB $(CUDAINCLUDE)
LDFLAGS   := 

OBJS = \
	sapporo.o \
	send_fetch_data.o \
	sapporoG6lib.o 

CUOBJS = host_evaluate_gravity.cu_o
TARGET = libsapporo.a

all:  $(CUDA_TK) $(TARGET) $(PROG)

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

