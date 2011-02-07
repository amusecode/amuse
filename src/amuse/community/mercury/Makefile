
MPIF90 ?= mpif90

FC = $(MPIF90)
SRC=src

AMUSE_DIR?=../../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py

LIBS=$(SRC)/libmercury.a

PIC=#-fPIC

FFLAGSOLD = -O -ip -u -g -traceback $(PIC)

FFLAGS2= -I./src/

CLASSNAME=Mercury

GLCLASSNAME=glMercury

all: mercury_worker

%.o: %.f90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2)  -c -o $@ $< 

%.o: %.F90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) -c -o $@ $< 
        
mercury_worker: worker.f90 interface.o
	$(FC) $(FFLAGS) $(FFLAGS2) $^ -o $@ $(LIBS)

mercury_worker_gl: glworker.f90 interface.o
	make -C src gl FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(FC) $(FFLAGS) $^ -o $@  $(GLOBJECTS) \
	$(GLLIB)  $(X11LIB) $(THREADLIB) $(LIBS)

worker.f90: interface.py
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

glworker.f90: interface.py
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(CODE_GENERATOR) --type=f90 $< $(GLCLASSNAME) -o $@

clean:
	make -C src clean
	rm -f *.o *.pyc *.bck worker worker.f90 glworker glworker.f90
	rm -f mercury_worker mercury_worker_gl
