CXX ?= g++
CC  ?= gcc

.SUFFIXES: .o .f90 .cpp .cu

CUDA_SDK ?= /disks/koppoel1/CUDA23/cudasdk/C
CUDA_TK  ?= /disks/koppoel1/CUDA23/cuda

OFLAGS = -O0 -g -Wall 
OFLAGS = -O3 -g -mpreferred-stack-boundary=4 -funroll-loops -fforce-addr
CFLAGS ?= $(OFLAGS)  -I$(CUDA_TK)/include
CXXFLAGS ?= $(CFLAGS)

NVCC      ?= $(CUDA_TK)/bin/nvcc

NVCCFLAGS ?= -D_DEBUG -m64 --compiler-options -fno-inline \
		--maxrregcount=32 -O0 -g -I$(CUDA_SDK)/common/inc

#device emulation version:
#NVCCFLAGS = -D_DEBUG --device-emulation --compiler-options -fno-inline \
#		--maxrregcount=32 -O0 -g -I$(CUDA_SDK)/common/inc

#muse version:
#NVCCFLAGS = -D_DEBUG --device-emulation \
#	 -Xcompiler "-fPIC" --maxrregcount=32  -O3 -I$(CUDA_SDK)/common/inc


#PROG = nbint
PROG = test_gravity

OCTLIB = octgrav

LIBOBJS = load_bodies.o \
	evaluate_gravity.o \
	reorder_bodies.o \
	memory.o \
	device_data.o \
	CUDA_evaluate_gravity.o

OBJS = $(PROG).o 

CUOBJS = host_evaluate_gravity.cu_o

CUDA_LIBDIRS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64 
CUDA_LIBS ?= -lcudart

LIBS = $(CUDA_LIBDIRS) $(CUDA_LIBS)

all: $(OCTLIB) $(PROG) 

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@ -L. -l$(OCTLIB)

$(OCTLIB): $(LIBOBJS) $(CUOBJS)
	ar qv lib$@.a $^
	ranlib lib$@.a

.cpp.o: 
	$(CXX) $(CXXFLAGS) -c $< -o $@

.f90.o:
	$(F90) $(F90FLAGS) -c $< -o $@

%.cu_o:  %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

convert:
	ifort -O3 -o convert convert.f90

clean:
	/bin/rm -rf *.o *.cu_o
	/bin/rm -rf $(PROG) lib$(OCTLIB).a

$(PROG).o:  octgrav.h
load_data.o:  octgrav.h
build_balanced_tree.o: octgrav.h
compute_multipole_moments.o: octgrav.h
sort_leaves.o: octgrav.h
reorder_bodies.o: octgrav.h
evaluate_gravity.o: octgrav.h
memory.o: octgrav.h
prepare_data_for_device.o: octgrav.h
copy.o: octgrav.h
CUDA_evaluate_gravity.o: octgrav.h

host_evaluate_gravity.cu_o: dev_evaluate_gravity.cu dev_octgrav_tex.cuh 
