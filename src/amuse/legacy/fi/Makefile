FC = mpif90

SRC=src

LIBS=$(SRC)/libfi.a

GLOBJECTS=$(SRC)/fthread.o $(SRC)/stars3d.o $(SRC)/libpt.a

GLINCLUDE = -I/usr/include/GL -I/home/inti/libraries/f90gl-1.2.11/include/GL/

GLLIB = -L/home/inti/libraries/f90gl-1.2.11/lib \
-lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut  -lguide 

X11LIB = -L/usr/X11R6/lib -lXaw -lXt -lXmu -lXi -lXext -lX11

FFTWLIB = -lfftw3 -lfftw3_threads

THREADLIB = -L. -lpthread src/libpt.a

OMPFLAGS= # -openmp -openmp_report0 
PIC=#-fPIC

# gfortran flags

ifeq ($(FORTRAN), gfortran)
FFLAGS = -fdefault-real-8 -O3 -ffree-line-length-256 -frecord-marker=4 $(PIC)
endif

ifeq  ($(FORTRAN), ifort)
# ifort flags
FFLAGS= -r8 -O -ip -u $(PIC)
endif

#FFLAGS= -fdefault-real-8 -O -ffree-line-length-256 -frecord-marker=4 $(PIC)


FFLAGS+= $(OMPFLAGS) 

FFLAGS+= -I./src/

LIBS+=$(FFTWLIB)

all: worker

%.o: %.f90 Makefile
	$(FC) $(FFLAGS)  -c -o $@ $< 
        
worker: worker.f90 interface.o
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS)

glworker: glworker.f90 interface.o
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"
	make -C src gl FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(FC) $(FFLAGS) $^ -o $@  $(GLOBJECTS) \
	$(GLLIB)  $(X11LIB) $(THREADLIB) $(LIBS)

worker.f90: interface.py __init__.py
	../../../../bin/create_f90_worker.py $< fi > $@

glworker.f90: interface.py __init__.py
	../../../../bin/create_f90_worker.py $< glfi > $@

__init__.py:
	touch $@

clean:
	make -C src purge
	rm -f __init__.py
	rm -f *.o *.pyc *.bck worker worker.f90 glworker glworker.f90
