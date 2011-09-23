MPIF90 ?= mpif90
FORTRAN ?= gfortran
FC ?= gfortran

SRC=src

AMUSE_DIR?=../../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py

LIBS=$(SRC)/libsphray.a

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi

FSLIBS ?= -L$(AMUSE_DIR)/lib/forsockets -lforsockets
FSINC ?= -I$(AMUSE_DIR)/lib/forsockets

PIC=#-fPIC

ifeq ($(findstring gfortran, $(notdir $(FORTRAN))), gfortran)
FFLAGS= -cpp -O5
endif

ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
IFORTFLAGS= -fpe0 -g -traceback -heap-arrays -fpp -u
FFLAGS= $(IFORTFLAGS) -O3 -no-prec-div -p
endif

FFLAGS2= -I./src/

CLASSNAME=SPHRayInterface

all: sphray_worker

%.o: %.f90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2)  -c -o $@ $< 

%.o: %.F90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) -c -o $@ $< 

sphray_worker: sphray_worker.f90 interface.o
	$(MPIF90) $(FFLAGS) $(FFLAGS2) $^ -o $@ $(LIBS) $(SCLIBS)
	
sphray_worker_sockets: sphray_worker-sockets.f90 interface.o
	$(FC) $(FFLAGS) $(FFLAGS2) $(FSINC) $^ -o $@ $(LIBS) $(SCLIBS) $(FSLIBS)

sphray_worker.f90: interface.py
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@
	
sphray_worker-sockets.f90: interface.py
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(CODE_GENERATOR) --type=f90 --mode=sockets $< $(CLASSNAME) -o $@

clean:
	make -C src clean
	rm -f *.o *.pyc *.bck
	rm -f sphray_worker sphray_worker.f90
	rm -f sphray_worker_sockets sphray_worker-sockets.f90
