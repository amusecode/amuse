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

LIBS+=$(SRC)/libfractal.a

PIC=#-fPIC

ifeq ($(findstring gfortran, $(notdir $(FORTRAN))), gfortran)
FFLAGS= $(FCFLAGS) -fdefault-real-8 -cpp -O3
endif

ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
IFORTFLAGS= -r8 -u
FFLAGS= $(FCFLAGS) $(IFORTFLAGS) -O2 
endif

FFLAGS2= -I./src/

CLASSNAME=FractalClusterInterface

all: fractal_worker

%.o: %.f90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2)  -c -o $@ $<

%.o: %.F90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) -c -o $@ $<

fractal_worker: fractal_worker.f90 interface.o
	$(MPIF90) $(FFLAGS) $(FFLAGS2) $(FS_FLAGS) $(LDFLAGS) $^ -o $@ $(LIBS) $(FS_LIBS)  $(LIBS)

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

