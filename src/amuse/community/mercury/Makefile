# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPIF90 ?= openmpif90

FC = $(MPIF90)
SRC=src

CODE_GENERATOR = $(AMUSE_DIR)/build.py

LIBS=$(SRC)/libmercury.a

SCINC ?= -I$(AMUSE_DIR)/lib/stopcond
SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi

FSLIBS ?= -L$(AMUSE_DIR)/lib/forsockets -lforsockets
FSINC ?= -I$(AMUSE_DIR)/lib/forsockets

PIC=#-fPIC

FFLAGSOLD = -O -ip -u -g -traceback $(PIC)

FFLAGS2= -I./src/

CLASSNAME=MercuryInterface

GLCLASSNAME=glMercury

all: mercury_worker

%.o: %.f90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) $(SCINC)  -c -o $@ $< 

%.o: %.F90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) $(SCINC) -c -o $@ $< 

mercury_worker: worker.f90 interface.o
	$(FC) $(FFLAGS) $(FFLAGS2) $^ -o $@ $(LIBS) $(SCLIBS)
	
mercury_worker_sockets: worker-sockets.f90 interface.o
	$(FC) $(FFLAGS) $(FFLAGS2) $(FSINC) $^ -o $@ $(LIBS) $(SCLIBS) $(FSLIBS)

mercury_worker_gl: glworker.f90 interface.o
	make -C src gl FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(FC) $(FFLAGS) $^ -o $@  $(GLOBJECTS) \
	$(GLLIB)  $(X11LIB) $(THREADLIB) $(LIBS) $(SCLIBS)

worker.f90: interface.py
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@
	
worker-sockets.f90: interface.py
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(CODE_GENERATOR) --type=f90 --mode=sockets $< $(CLASSNAME) -o $@

glworker.f90: interface.py
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(CODE_GENERATOR) --type=f90 $< $(GLCLASSNAME) -o $@

clean:
	make -C src clean
	rm -f *.o *.pyc *.bck worker worker.f90 glworker glworker.f90
	rm -f mercury_worker mercury_worker_gl
	rm -f mercury_worker_sockets worker-sockets.f90
	rm -f forsockets.mod
	rm -f src/lib*.a
