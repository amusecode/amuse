
MPIF90 ?= mpif90
FORTRAN ?= gfortran
FC = $(MPIF90)

SRC=src

AMUSE_DIR?=../../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py

NOMPILIBS=$(SRC)/libcapreole.a
MPILIBS=$(SRC)/libcapreole_mpi.a

GLOBJECTS=$(SRC)/fthread.o $(SRC)/viewer.o $(SRC)/libpt.a

GLINCLUDE = -I/usr/include/GL -I/data/pelupes/libraries/f90gl-1.2.11/include/GL/

GLLIB = -L/data/pelupes/libraries/f90gl-1.2.11/lib \
-lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut  -lguide 

X11LIB = -L/usr/X11R6/lib -lXaw -lXt -lXmu -lXi -lXext -lX11

THREADLIB = -L. -lpthread src/libpt.a

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

OMPFLAGS= # -openmp -openmp_report0 
PIC=#-fPIC

ifeq ($(FORTRAN), gfortran)
FFLAGS += -fdefault-real-8 -O3 -ffree-line-length-256 -frecord-marker=4 $(PIC)
endif

ifeq  ($(FORTRAN), ifort)
# ifort flags
FFLAGS += -O -ip -u -g -traceback -fpp1 -DIFORT  $(PIC)
endif

MPI_FLAGS = -DMPI # 

FFLAGS+= $(OMPFLAGS) 

FFLAGS2= -I./src/

CLASSNAME=CapreoleInterface

GLCLASSNAME=GLCapreoleInterface

FSLIB = -L$(AMUSE_DIR)/lib/forsockets -lforsockets
FSINC = -I$(AMUSE_DIR)/lib/forsockets

all: capreole_worker

%.o: %.f90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) $(SCINC)  -c -o $@ $< 

%.o: %.F90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) $(SCINC) -c -o $@ $< 

$(MPILIBS):
	make -C src clean
	make -C src amuse_interface_mpi F90FLAGS1="$(FFLAGS)" FC="$(FC)"
	
$(NOMPILIBS):
	make -C src clean
	make -C src amuse_interface F90FLAGS1="$(FFLAGS)" FC="$(FC)"
	
capreole_worker: $(NOMPILIBS) worker.f90 interface.o
	$(FC) $(FFLAGS) $(FFLAGS2) $(SCINC) $^ -o $@ $(NOMPILIBS) $(SCLIBS)
	
capreole_worker_sockets: $(NOMPILIBS) worker-sockets.F90 interface.o
	$(FC) $(FFLAGS) $(FFLAGS2) $(SCINC) $(FSINC) $^ -o $@ $(NOMPILIBS) $(SCLIBS) $(FSLIB)

capreole_worker_mpi: FFLAGS2+=$(MPI_FLAGS)
capreole_worker_mpi: $(MPILIBS) worker_mpi.f90 interface.o
	$(FC) $(FFLAGS) $(FFLAGS2) $(SCINC) $^ -o $@ $(MPILIBS) $(SCLIBS)

capreole_worker_gl: $(NOMPILIBS) glworker.f90 interface.o
	make -C src clean
	make -C src gl FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(FC) $(FFLAGS) $(SCINC) $^ -o $@  $(GLOBJECTS) \
	$(GLLIB)  $(X11LIB) $(THREADLIB) $(NOMPILIBS) $(SCLIBS)

worker.f90: interface.py 
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@
	
worker-sockets.F90: interface.py 
	$(CODE_GENERATOR) --type=f90 --mode=sockets $< $(CLASSNAME) -o $@

worker_mpi.f90: interface.py 
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

glworker.f90: interface.py 
	make -C src amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)"
	$(CODE_GENERATOR) --type=f90 $< $(GLCLASSNAME) -o $@

clean:
	make -C src clean
	rm -f *.o *.pyc *.bck worker worker.f90 glworker glworker.f90 worker-sockets.F90 worker-sockets.f90
	rm -f capreole_worker capreole_worker_mpi capreole_worker_gl worker_mpi.f90 capreole_worker_sockets
