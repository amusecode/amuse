# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../..
-include ${AMUSE_DIR}/config.mk


CXX ?= g++
CC  ?= gcc
LD  ?= g++ 

.SUFFIXES: .o .cpp .ptx .cu

CUDA_TK  ?= /usr/local/cuda

# OFLAGS = -O0 -g -Wall 
#OFLAGS = -O3 -g -Wall -fopenmp  -D_OCL_
#OFLAGS = -O3 -g -Wall -fopenmp
#OFLAGS = -g -Wall -fopenmp -D_OCL_
OFLAGS = -g -Wall -fopenmp 

CXXFLAGS =  -fPIC $(OFLAGS) -I$(CUDA_TK)/include  $(CL_CFLAGS)


# NVCC      = $(CUDA_TK)/bin/nvcc  --device-emulation
# NVCCFLAGS = -D_DEBUG -O0 -g -I$(CUDA_SDK)/common/inc -arch=sm_12 --maxrregcount=64  --opencc-options -OPT:Olimit=0 -I$(CUDPP)/cudpp/include
NVCC      = $(CUDA_TK)/bin/nvcc  

#NVCCFLAGS = -arch sm_12  #<-- gives slightly faster kernels, because of limited DP support
#NVCCFLAGS ?= -arch sm_20
NVCCFLAGS ?= $(NVCC_FLAGS) -arch=sm_12

# Use with Mac OS X
# NVCCFLAGS = -arch sm_12 -Xcompiler="-Duint=unsigned\ int"

LDFLAGS = -lcuda $(CL_LIBS) -fopenmp 


INCLUDEPATH = ./include
CXXFLAGS  += -I$(INCLUDEPATH) -D__INCLUDE_KERNELS__
NVCCFLAGS += -I$(INCLUDEPATH)

CUDAKERNELSPATH = CUDA
CUDAKERNELS = kernels4th.cu  kernels4thDP.cu kernels6thDP.cu kernelsG5DS.cu kernelsG5SP.cu

CUDAPTX = $(CUDAKERNELS:%.cu=$(CUDAKERNELSPATH)/%.ptx)

SRCPATH = src
SRC = sapporohostclass.cpp sapporoG6lib.cpp sapporoYeblib.cpp sapporoG5lib.cpp sapporo6thlib.cpp
OBJ = $(SRC:%.cpp=%.o)

PROG = 
#PROG = main

#LIBOBJ = sapporohostclass.o sapporoG6lib.o sapporoYeblib.o sapporoG5lib.o
LIBOBJ = sapporohostclass.o sapporoG6lib.o sapporoG5lib.o sapporo6thlib.o
TARGET = libsapporo.a


all:	  $(OBJ) $(CUDAPTX) $(PROG) $(TARGET)
#all: 	  $(OBJ) $(PROG)
kernels:  $(CUDAPTX)

$(PROG): $(OBJ)
	$(LD) $(LDFLAGS) $^ -o $@ 

$(TARGET): $(LIBOBJ)
	ar qv $@ $^        

%.o: $(SRCPATH)/%.cpp
	$(CXX) $(CXXFLAGS)  -c $< -o $@

$(CUDAKERNELSPATH)/%.ptx: $(CUDAKERNELSPATH)/%.cu
	$(NVCC) $(NVCCFLAGS) -ptx $< -o $@

clean:
	/bin/rm -rf *.o *.ptx *.a CUDA/*ptx
	rm *.ptxh include/*.ptxh

$(OBJ): $(INCLUDEPATH)/*.h

#testkernel.ptx: 

ifneq (,findstring($(CUDA_ENABLED),"yes"))
sapporohostclass.o : include/kernels4th.ptxh include/kernels6thDP.ptxh include/kernelsG5DS.ptxh include/kernelsG5SP.ptxh sapporohostclass.h sapdevclass.h $(INCLUDEPATH)/defines.h
$(CUDAKERNELSPATH)/kernels4th.ptx : $(INCLUDEPATH)/defines.h
$(CUDAKERNELSPATH)/kernels4thDP.ptx : $(INCLUDEPATH)/defines.h
$(CUDAKERNELSPATH)/kernels6thDP.ptx : $(INCLUDEPATH)/defines.h
$(CUDAKERNELSPATH)/kernelsG5DS.ptx : $(INCLUDEPATH)/defines.h
$(CUDAKERNELSPATH)/kernelsG5SP.ptx : $(INCLUDEPATH)/defines.h
endif


libsapporo.a : sapporohostclass.o

include/%.ptxh:CUDA/%.ptx
	xxd -i $< $@






