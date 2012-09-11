# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPIFC ?= mpif90
FORTRAN ?= gfortran
FC = $(MPIFC)

SRC=src
BUILDDIR = build
MPI_BUILDDIR = build_mpi

CODE_GENERATOR = $(AMUSE_DIR)/build.py

NOMPILIBS=$(BUILDDIR)/libcapreole.a
MPILIBS=$(MPI_BUILDDIR)/libcapreole_mpi.a

GLOBJECTS=$(SRC)/fthread.o $(SRC)/viewer.o $(SRC)/libpt.a

GLINCLUDE = -I/usr/include/GL -I/home/inti/libraries/f90gl-1.2.11/include/GL/

GLLIB = -L/home/inti/libraries/f90gl-1.2.11/lib \
-lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut  

X11LIB = -L/usr/X11R6/lib64 -lXaw -lXt -lXmu -lXi -lXext -lX11

THREADLIB = -L. -lpthread src/libpt.a

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcond
SCMPILIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi
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

FFLAGS2= 

CLASSNAME=CapreoleInterface

GLCLASSNAME=GLCapreoleInterface

FSLIB = -L$(AMUSE_DIR)/lib/forsockets -lforsockets
FSINC = -I$(AMUSE_DIR)/lib/forsockets

all: capreole_worker

$(BUILDDIR)/%.o: %.F90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) -I$(BUILDDIR) $(SCINC) -c -o $@ $< 
	

$(MPI_BUILDDIR)/%.o: %.F90 Makefile
	$(FC) $(FFLAGS) $(FFLAGS2) $(MPI_FLAGS) -I$(MPI_BUILDDIR) $(SCINC) -c -o $@ $< 

$(BUILDDIR)/Makefile:
	-mkdir $(BUILDDIR)
	cp $(SRC)/Makefile $(BUILDDIR)/Makefile

$(NOMPILIBS): $(BUILDDIR)/Makefile
	make -C $(BUILDDIR) amuse_interface  VPATH=../src F90FLAGS1="$(FFLAGS)" FC="$(FC)"
	
$(MPI_BUILDDIR)/Makefile:
	-mkdir $(MPI_BUILDDIR)
	cp $(SRC)/Makefile $(MPI_BUILDDIR)/Makefile
	
$(MPILIBS): $(MPI_BUILDDIR)/Makefile
	make -C $(MPI_BUILDDIR) amuse_interface_mpi  VPATH=../src F90FLAGS1="$(FFLAGS)" FC="$(FC)" MPIFC="$(MPIFC)"
	
capreole_worker: $(NOMPILIBS) worker_code.f90 $(BUILDDIR)/interface.o
	$(FC) $(FFLAGS) $(FFLAGS2)  -I$(BUILDDIR) $(SCINC) $^ -o $@ $(NOMPILIBS) $(SCLIBS)
	
capreole_worker_sockets: $(NOMPILIBS) worker_code_sockets.F90 $(BUILDDIR)/interface.o
	$(FC) $(FFLAGS) $(FFLAGS2) $(SCINC) $(FSINC) $^ -o $@ $(NOMPILIBS) $(SCLIBS) $(FSLIB)

capreole_worker_mpi: $(MPILIBS) worker_code_mpi.f90 $(MPI_BUILDDIR)/interface.o
	$(FC) $(FFLAGS) $(FFLAGS2) $(MPI_FLAGS) -I$(MPI_BUILDDIR) $(SCINC) $^ -o $@ $(MPILIBS) $(SCMPILIBS)

#capreole_worker_gl: $(NOMPILIBS) glworker.f90 $(BUILDDIR)/interface.o
#	make -C src clean
#	make -C src gl FFLAGS="$(FFLAGS)" FC="$(FC)"
#	$(FC) $(FFLAGS) $(SCINC) $^ -o $@  $(GLOBJECTS) \
#	$(GLLIB)  $(X11LIB) $(THREADLIB) $(NOMPILIBS) $(SCLIBS)

worker_code.f90: interface.py 
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@
	
worker_code_sockets.F90: interface.py 
	$(CODE_GENERATOR) --type=f90 --mode=sockets $< $(CLASSNAME) -o $@

worker_code_mpi.f90: interface.py 
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

clean:
	make -C src clean
	$(RM) -f *.o *.pyc *.bck worker worker.f90 glworker glworker.f90 worker-sockets.F90 worker-sockets.f90
	$(RM) -f worker_code.f90 worker_code_sockets.F90 worker_code_mpi.f90
	$(RM) -f capreole_worker capreole_worker_mpi capreole_worker_gl worker_mpi.f90 capreole_worker_sockets
	$(RM) -Rf build build_mpi
	$(RM) -f forsockets.mod
