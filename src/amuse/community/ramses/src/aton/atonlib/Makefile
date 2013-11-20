# Makefile for libaton

# Grid size
# ---------
#
# The ATON grid size is fixed at compile time.
# You can set the dimensions here.
# For small tests, use 32x32x32.
# For a 256^3 simulation on 16 nodes, use 128x128x64.
NCELLX=32
NCELLY=32
NCELLZ=32
# TODO: Allow the grid size to be set at runtime.

# CUDA library paths
# ------------------
#
# These may need to be changed for your setup.
# For Titane:
LIBS = -lm -lcuda  -L/opt/cuda/lib -lcudart -I/opt/cuda/include
# For Mayall:
#LIBS = -lm -lcuda  -L/usr/local/cuda/lib64 -L/home/itp/tstranex/cuda/sdk/C/lib -lcudart -I/usr/local/cuda/include

# ---------------------------------------

CC_OBJS= cpu_numerical.o cpu_memory.o validation.o aton_fortran.o gpu_init.o aton_cpp.o
CUDA_OBJS= gpu_memory.o gpu_numerical.o
OBJECTS = $(CUDA_OBJS) $(CC_OBJS)
VPATH=internal

DEFINES = -DS_100000 -DNCELLX=${NCELLX} -DNCELLY=${NCELLY} -DNCELLZ=${NCELLZ} # -DTEST7_RAYTRACING_TRICKS
PUBLIC_INCLUDES=$(shell pwd)


NVCC= nvcc -arch=sm_13 #--device-emulation -g
CC = mpic++ -O2
F90 = mpif90
# For todi.cscs.ch, load module PrgEnv-gnu and set:
# CC = cc -O2
# F90 = ftn


EXECUTABLE = 

libaton.a: $(OBJECTS)
	ar rcs libaton.a $(OBJECTS)

aton_cpu_test: $(OBJECTS) internal/aton_cpu_test.o
	$(F90) -o aton_cpu_test internal/aton_cpu_test.o $(OBJECTS) $(LIBS)

stromgren_test: $(OBJECTS) examples/stromgren_test.o
	$(CC) -o stromgren_test stromgren_test.o $(OBJECTS) $(LIBS)

%.o:%.f90
	$(F90) -c $*.f90
.cc.o:
	$(CC) $(DEFINES) $(LIBS) -I$(PUBLIC_INCLUDES) -c $^
%.o:%.cu
	$(NVCC) $(DEFINES) $(LIBS) -I. -c $^

.PHONY : all
all: $(CUDA_OBJS)

clean:
	rm -f *.o *.cudafe1.* *.cudafe2.* *.hash *.ptx *fatbin.c *.cubin *.cpp* $(EXECUTABLE) libaton.a stromgren_test *~ *.mod
