FC = mpif90

SRC=src

AMUSE_DIR?=../../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py

NOMPILIBS=$(SRC)/libcapreole.a
MPILIBS=$(SRC)/libcapreole_mpi.a

GLOBJECTS=$(SRC)/fthread.o $(SRC)/viewer.o $(SRC)/libpt.a

GLINCLUDE = -I/usr/include/GL -I/home/inti/libraries/f90gl-1.2.11/include/GL/

GLLIB = -L/home/inti/libraries/f90gl-1.2.11/lib \
-lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut  -lguide 

X11LIB = -L/usr/X11R6/lib -lXaw -lXt -lXmu -lXi -lXext -lX11

THREADLIB = -L. -lpthread src/libpt.a

OMPFLAGS= # -openmp -openmp_report0 
PIC=#-fPIC

ifeq ($(FORTRAN), gfortran)
FFLAGS = -fdefault-real-8 -O3 -ffree-line-length-256 -frecord-marker=4 $(PIC)
endif

ifeq  ($(FORTRAN), ifort)
# ifort flags
FFLAGS= -O -ip -u -g -traceback -DIFORT  $(PIC)
endif

MPI_FLAGS = -DMPI # 

FFLAGS+= $(OMPFLAGS) 

FFLAGS2= -I./src/

CLASSNAME=capreole

all: worker

%.o: %.f90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2)  -c -o $@ $< 

%.o: %.F90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) -c -o $@ $< 
        
worker: worker.f90 interface.o
	$(FC) $(FFLAGS) $(FFLAGS2) $^ -o $@ $(NOMPILIBS)

worker_mpi: FFLAGS2+=$(MPI_FLAGS)
worker_mpi: worker_mpi.f90 interface.o
	$(FC) $(FFLAGS) $(FFLAGS2) $^ -o $@ $(MPILIBS)

glworker: glworker.f90 interface.o
	make -C src gl FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(FC) $(FFLAGS) $^ -o $@  $(GLOBJECTS) \
	$(GLLIB)  $(X11LIB) $(THREADLIB) $(NOMPILIBS)

worker.f90: interface.py __init__.py
	make -C src amuse_interface F90FLAGS1="$(FFLAGS)" FC="$(FC)"
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@


worker_mpi.f90: interface.py __init__.py
	make -C src amuse_interface_mpi F90FLAGS1="$(FFLAGS)" FC="$(FC)"
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

glworker.f90: interface.py __init__.py
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

__init__.py:
	touch $@

clean:
	make -C src clean
	rm -f __init__.py
	rm -f *.o *.pyc *.bck worker worker.f90 glworker glworker.f90
	rm -f worker_mpi worker_mpi.f90
