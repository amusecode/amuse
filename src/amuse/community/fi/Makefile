
MPIF90 ?= mpif90
FORTRAN ?= gfortran
FC = $(MPIF90)
GL_PATH ?= /home/inti/libraries/f90gl-1.2.11

SRC=src
BUILDDIR = build
PERIODIC_BUILDDIR = build_periodic

AMUSE_DIR?=../../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py



GLOBJECTS=$(BUILDDIR)/fthread.o $(BUILDDIR)/stars3d.o $(BUILDDIR)/libpt.a
GLINCLUDE = -I/usr/include/GL -I$(GL_PATH)/include/GL/
GLLIB = -L$(GL_PATH)/lib -lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut

X11LIB = -L/usr/X11R6/lib -lXaw -lXt -lXmu -lXi -lXext -lX11

FFTW_LIBS ?= -lfftw3 

THREADLIB = -L. -lpthread $(BUILDDIR)/libpt.a

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

FSLIBS ?= -L$(AMUSE_DIR)/lib/forsockets -lforsockets
FSINC ?= -I$(AMUSE_DIR)/lib/forsockets

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
GLLIB = $(GLLIB) -lguide 
endif

#FFLAGS= -fdefault-real-8 -O -ffree-line-length-256 -frecord-marker=4 $(PIC)


FFLAGS+= $(OMPFLAGS) 

FFLAGS+= -I.

LIBS+=$(FFTW_LIBS)
LIBS+=-lfftw3_threads

FILES=$(SRC)/muse_helpers.f90 $(SRC)/mem.f90 $(SRC)/elements.f90 $(SRC)/cooling.f90 \
  $(SRC)/ionize.f90 $(SRC)/molecules.f90 $(SRC)/H2cool.f90 $(SRC)/starprop.f90 \
  $(SRC)/random.f90 $(SRC)/io.f90 $(SRC)/io-old.f $(SRC)/diagnostics.f90 $(SRC)/readparam.f90 \
  $(SRC)/ppread.f90 $(SRC)/blackhole.f90 $(SRC)/stepsystem.f90 $(SRC)/buildnearlist.f90 \
  $(SRC)/pmgrav.f90 $(SRC)/search.f90 $(SRC)/feedback.f90 $(SRC)/fixedhalo.f90 $(SRC)/init.f90 \
  $(SRC)/timestep.f90 $(SRC)/buildtree.f90 $(SRC)/density.f90 $(SRC)/io2.f90 $(SRC)/heating.f90 \
  $(SRC)/entdot.f90 $(SRC)/heco.f90 $(SRC)/ethstep.f90 $(SRC)/stellar.f90 $(SRC)/fcco.f90 \
  $(SRC)/project.f90 $(SRC)/clean.f90 $(SRC)/sort.f90 $(SRC)/util.f90 $(SRC)/gravity.f90 \
  $(SRC)/tidalsum.f90 $(SRC)/treewalk.f90 $(SRC)/gravsum.f90 $(SRC)/fuvflux.f90

AMUSE_INTERFACE_DEPS=../interface/gd.py ../interface/common.py ../../support/codes/core.py ../../support/interface.py


all: fi_worker

%.o: %.f90 Makefile
	$(FC) $(FFLAGS)  -c -o $@ $< 
	
$(BUILDDIR)/Makefile: 
	-mkdir $(BUILDDIR)
	ln -s ../src/Makefile  $(BUILDDIR)/Makefile
	
$(BUILDDIR)/libfi.a: $(BUILDDIR)/Makefile $(FILES) Makefile
	-mkdir $(BUILDDIR)
	make -C $(BUILDDIR) amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)" VPATH=../src 

fi_worker: worker.f90 interface.o $(BUILDDIR)/libfi.a
	$(FC) $(FFLAGS) $(SCINC) worker.f90 interface.o -o $@ $(BUILDDIR)/libfi.a  $(LIBS) $(SCLIBS)

fi_worker_sockets: worker-sockets.f90 interface.o $(BUILDDIR)/libfi.a
	$(FC) $(FFLAGS) $(SCINC) $(FSINC) worker-sockets.f90 interface.o -o $@ $(BUILDDIR)/libfi.a  $(LIBS) $(SCLIBS) $(FSLIBS)

fi_worker_gl: glworker.f90 interface.o $(BUILDDIR)/libfi.a
	make -C $(BUILDDIR) gl FFLAGS="$(FFLAGS)" FC="$(FC)" GLINCLUDE="$(GLINCLUDE)" GLLIB="$(GLLIB)" VPATH=../src 
	$(FC) $(FFLAGS) $(SCINC) glworker.f90 interface.o -o $@  $(GLOBJECTS) \
	$(GLLIB)  $(X11LIB) $(THREADLIB) $(BUILDDIR)/libfi.a $(LIBS) $(SCLIBS)

$(PERIODIC_BUILDDIR)/Makefile:
	-mkdir $(PERIODIC_BUILDDIR)
	ln -s ../src/Makefile $(PERIODIC_BUILDDIR)/Makefile
	
$(PERIODIC_BUILDDIR)/libfi.a: $(PERIODIC_BUILDDIR)/Makefile 
	./patch_globals_dot_h.py --periodic
	make -C $(PERIODIC_BUILDDIR) amuse_interface FFLAGS="-I./$(PERIODIC_BUILDDIR) $(FFLAGS)" FC="$(FC)" VPATH=../src 
	./patch_globals_dot_h.py
	
fi_worker_periodic: worker.f90 interface.po $(PERIODIC_BUILDDIR)/libfi.a
	$(FC) -I./$(PERIODIC_BUILDDIR) $(FFLAGS) $(SCINC) worker.f90 interface.po -o $@  $(PERIODIC_BUILDDIR)/libfi.a $(SCLIBS) $(FFTW_LIBS)

worker.f90: interface.py $(AMUSE_INTERFACE_DEPS)
	$(CODE_GENERATOR) --type=f90 $< FiInterface -o $@

worker-sockets.f90: interface.py $(AMUSE_INTERFACE_DEPS)
	$(CODE_GENERATOR) --type=f90 --mode=sockets $< FiInterface -o $@

glworker.f90: interface.py $(AMUSE_INTERFACE_DEPS)
	$(CODE_GENERATOR) --type=f90 $< GlFiInterface -o $@

interface.po:interface.f90
	$(FC) $(FFLAGS) -I./$(PERIODIC_BUILDDIR) -c -o $@ $<
    
clean:
	make -C src purge
	rm -f *.o *.pyc *.bck worker worker.f90 glworker glworker.f90
	rm -f fi_worker_periodic fi_worker_gl fi_worker
	rm -f fi_worker_sockets worker-sockets.f90	
	rm -f worker_code_periodic worker_code
	rm -Rf $(BUILDDIR)
	rm -Rf $(PERIODIC_BUILDDIR)

