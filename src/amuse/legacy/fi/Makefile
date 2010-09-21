
MPIF90 ?= mpif90
FORTRAN ?= gfortran
FC = $(MPIF90)

SRC=src

AMUSE_DIR?=../../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py

LIBS=$(SRC)/libfi.a

GLOBJECTS=$(SRC)/fthread.o $(SRC)/stars3d.o $(SRC)/libpt.a

GLINCLUDE = -I/usr/include/GL -I/home/inti/libraries/f90gl-1.2.11/include/GL/

GLLIB = -L/home/inti/libraries/f90gl-1.2.11/lib \
-lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut  -lguide 

X11LIB = -L/usr/X11R6/lib -lXaw -lXt -lXmu -lXi -lXext -lX11

FFTWLIB = -lfftw3 -lfftw3_threads

THREADLIB = -L. -lpthread src/libpt.a

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

OMPFLAGS=# -openmp -openmp_report0 
PIC=#-fPIC

# gfortran flags


ifeq ($(findstring gfortran, $(notdir $(FORTRAN))), gfortran)
FFLAGS = -fdefault-real-8 -O3 -ffree-line-length-256 -frecord-marker=4 $(PIC)
endif

ifeq ($(findstring f95, $(notdir $(FORTRAN))), f95)
FFLAGS = -fdefault-real-8 -O3 -ffree-line-length-256 -frecord-marker=4 $(PIC)
endif

ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
# ifort flags
FFLAGS= -r8 -O -ip -u $(PIC)
endif

#FFLAGS= -fdefault-real-8 -O -ffree-line-length-256 -frecord-marker=4 $(PIC)


FFLAGS+= $(OMPFLAGS) 

FFLAGS+= -I./src/

LIBS+=$(FFTWLIB)

FILES=$(SRC)/muse_helpers.f90 $(SRC)/mem.f90 $(SRC)/elements.f90 $(SRC)/cooling.f90 \
  $(SRC)/ionize.f90 $(SRC)/molecules.f90 $(SRC)/H2cool.f90 $(SRC)/starprop.f90 \
  $(SRC)/random.f90 $(SRC)/io.f90 $(SRC)/io-old.f $(SRC)/diagnostics.f90 $(SRC)/readparam.f90 \
  $(SRC)/ppread.f90 $(SRC)/blackhole.f90 $(SRC)/stepsystem.f90 $(SRC)/buildnearlist.f90 \
  $(SRC)/pmgrav.f90 $(SRC)/search.f90 $(SRC)/feedback.f90 $(SRC)/fixedhalo.f90 $(SRC)/init.f90 \
  $(SRC)/timestep.f90 $(SRC)/buildtree.f90 $(SRC)/density.f90 $(SRC)/io2.f90 $(SRC)/heating.f90 \
  $(SRC)/entdot.f90 $(SRC)/heco.f90 $(SRC)/ethstep.f90 $(SRC)/stellar.f90 $(SRC)/fcco.f90 \
  $(SRC)/project.f90 $(SRC)/clean.f90 $(SRC)/sort.f90 $(SRC)/util.f90 $(SRC)/gravity.f90 \
  $(SRC)/tidalsum.f90 $(SRC)/treewalk.f90 $(SRC)/gravsum.f90 $(SRC)/fuvflux.f90

AMUSE_INTERFACE_DEPS=../interface/gd.py ../interface/common.py ../../support/legacy/core.py ../../support/interface.py


all: worker

%.o: %.f90 Makefile
	$(FC) $(FFLAGS)  -c -o $@ $< 

$(SRC)/libfi.a: $(FILES) Makefile
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"

worker: worker.f90 interface.o $(SRC)/libfi.a
	$(FC) $(FFLAGS) $(SCINC) worker.f90 interface.o -o $@ $(LIBS) $(SCLIBS)

glworker: glworker.f90 interface.o $(SRC)/libfi.a
	make -C src gl FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(FC) $(FFLAGS) $(SCINC) glworker.f90 interface.o -o $@  $(GLOBJECTS) \
	$(GLLIB)  $(X11LIB) $(THREADLIB) $(LIBS) $(SCLIBS)

worker.f90: interface.py __init__.py $(AMUSE_INTERFACE_DEPS)
	$(CODE_GENERATOR) --type=f90 $< FiInterface -o $@

glworker.f90: interface.py __init__.py $(AMUSE_INTERFACE_DEPS)
	$(CODE_GENERATOR) --type=f90 $< GlFiInterface -o $@


__init__.py:
	touch $@

clean:
	make -C src purge
	rm -f __init__.py
	rm -f *.o *.pyc *.bck worker worker.f90 glworker glworker.f90
