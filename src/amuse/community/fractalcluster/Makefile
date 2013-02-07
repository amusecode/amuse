# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPIF90 ?= mpif90
FORTRAN ?= gfortran
FC ?= gfortran

SRC=src

CODE_GENERATOR = $(AMUSE_DIR)/build.py

LIBS=$(SRC)/libfractal.a

PIC=#-fPIC

ifeq ($(findstring gfortran, $(notdir $(FORTRAN))), gfortran)
FFLAGS= -fdefault-real-8 -cpp -O3
endif

ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
IFORTFLAGS= -r8 -u
FFLAGS= $(IFORTFLAGS) -O2 
endif

FFLAGS2= -I./src/

CLASSNAME=FractalClusterInterface

all: fractal_worker

%.o: %.f90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2)  -c -o $@ $<

%.o: %.F90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) -c -o $@ $<

fractal_worker: fractal_worker.f90 interface.o
	$(MPIF90) $(FFLAGS) $(FFLAGS2) $(FS_FLAGS) $^ -o $@ $(LIBS) $(FS_LIBS)

fractal_worker.f90: interface.py
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

clean:
	make -C src clean
	$(RM) -f *.o *.pyc *.bck
	$(RM) -f fractal_worker fractal_worker.f90
	$(RM) -f fractal_worker_sockets fractal_worker-sockets.f90
	$(RM) -f forsockets.mod
	$(RM) -f amuse_worker_module.mod
	$(RM) -f src/lib*.a

