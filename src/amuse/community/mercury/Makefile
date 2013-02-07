# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPIFC ?= openmpif90

FC = $(MPIFC)
SRC=src

CODE_GENERATOR = $(AMUSE_DIR)/build.py

LIBS=$(SRC)/libmercury.a

PIC=#-fPIC

FFLAGSOLD = -O -ip -u -g -traceback $(PIC)

FFLAGS2= -I./src/

CLASSNAME=MercuryInterface

F77FLAGS = 
ifneq (,$(findstring xlf, $(notdir $(FORTRAN))))
F77FLAGS += -qfixed -qsuffix=f=for
endif


all: mercury_worker

%.o: %.f90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) $(SC_FLAGS)  -c -o $@ $< 

%.o: %.F90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) $(SC_FLAGS) -c -o $@ $< 

mercury_worker: worker.f90 interface.o
	$(FC) $(FFLAGS) $(FFLAGS2) $(SC_FLAGS) $(FS_FLAGS) $^ -o $@ $(LIBS) $(SC_FCLIBS) $(FS_LIBS)
	
worker.f90: interface.py
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)" F77FLAGS="$(F77FLAGS)"
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@
	
clean:
	make -C src clean
	rm -f *.o *.pyc *.bck worker worker.f90 glworker glworker.f90
	rm -f mercury_worker mercury_worker_gl
	rm -f mercury_worker_sockets worker-sockets.f90
	rm -f amuse_worker_module.mod
	rm -f amuseinterface.mod
	rm -f forsockets.mod
	rm -f src/lib*.a
