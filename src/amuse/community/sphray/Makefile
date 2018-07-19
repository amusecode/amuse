# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(strip $(AMUSE_DIR)),)
    AMUSE_DIR = ../../../..
endif
-include ${AMUSE_DIR}/config.mk

ifeq ($(strip $(MPIFC)),)
    MPIFC = mpif90
endif
ifeq ($(strip $(MPIF90)),)
    MPIF90 = $(MPIFC)
endif
ifeq ($(strip $(FORTRAN)),)
    FORTRAN = gfortran
endif
ifeq ($(strip $(FC)),)
    FC = gfortran
endif

SRC=src

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CODE_LIBS=$(SRC)/libsphray.a

PIC=#-fPIC

ifeq ($(findstring gfortran, $(notdir $(FORTRAN))), gfortran)
FFLAGS= $(FCFLAGS) -cpp -O0  # problems with ANY optimizations in gfortran 4.4 and below
endif

ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
IFORTFLAGS= -fpe0 -g -traceback -heap-arrays -fpp -u
FFLAGS= $(FCFLAGS) $(IFORTFLAGS) -O3 -no-prec-div -p
endif

FFLAGS2= -I./src/

CLASSNAME=SPHRayInterface

all: sphray_worker

%.o: %.f90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2)  -c -o $@ $< 

%.o: %.F90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) -c -o $@ $< 

sphray_worker: sphray_worker.f90 interface.o
	$(MPIF90) $(FFLAGS) $(FFLAGS2) $(FS_FLAGS) $(LDFLAGS) $^ -o $@ $(CODE_LIBS) $(SC_MPI_CLIBS) $(FS_LIBS) $(LIBS) 

sphray_worker.f90: interface.py
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

clean:
	make -C src clean
	$(RM) -f *.o *.pyc *.bck
	$(RM) -f sphray_worker sphray_worker.f90
	$(RM) -f sphray_worker_sockets sphray_worker-sockets.f90
	$(RM) -f forsockets.mod
	$(RM) -f amuse_worker_module.mod
	$(RM) -f src/lib*.a
