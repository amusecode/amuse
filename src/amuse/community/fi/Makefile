# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPIFC ?= mpif90
FORTRAN=$(FC)
PYTHON ?= python
GL_PATH ?= /home/inti/libraries/f90gl-1.2.11

SRC=src
BUILDDIR = build
PERIODIC_BUILDDIR = build_periodic
OPENMP_BUILDDIR = build_mp

CODE_GENERATOR = $(AMUSE_DIR)/build.py


GLOBJECTS=$(BUILDDIR)/fthread.o $(BUILDDIR)/fi_viewer.o $(BUILDDIR)/libpt.a
PGLOBJECTS=$(PERIODIC_BUILDDIR)/fthread.o $(PERIODIC_BUILDDIR)/fi_viewer.o $(PERIODIC_BUILDDIR)/libpt.a
GLINCLUDE = -I/usr/include/GL -I$(GL_PATH)/include/GL/
GLLIB = -L$(GL_PATH)/lib -lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut

X11LIB = -L/usr/lib64 -lXaw -lXt -lXmu -lXi -lXext -lX11

FFTW_LIBS ?= -lfftw3  -lfftw3_threads

THREADLIB = -L. -lpthread $(BUILDDIR)/libpt.a
PTHREADLIB = -L. -lpthread $(PERIODIC_BUILDDIR)/libpt.a

OPENMP_FCFLAGS ?= -openmp -openmp_report0 

PIC=#-fPIC

F77FLAGS =
# gfortran flags
$(info  $(shell echo $(notdir $(FORTRAN)) |cut -f1 -d- ))
$(info  $(findstring gfortran, $(shell echo $(notdir $(FORTRAN)) |cut -f1 -d- )))
ifneq (,$(findstring gfortran, $(notdir $(FORTRAN))))
FFLAGS = $(FCFLAGS) -fdefault-real-8 -O3 -ffree-line-length-256 -frecord-marker=4 $(PIC)
endif


ifeq ($(findstring f95, $(notdir $(FORTRAN))), f95)
FFLAGS = $(FCFLAGS) -fdefault-real-8  -O3 -ffree-line-length-256 -frecord-marker=4 $(PIC)
endif

ifeq ($(findstring xlf, $(notdir $(FORTRAN))), xlf)
FFLAGS = -qrealsize=8 -O2 $(PIC)
F77FLAGS += -qfixed
endif

ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
# ifort flags
FFLAGS= $(FCFLAGS) -r8 -O2 -ip -u $(PIC)
#GLLIB = $(GLLIB) -lguide 
endif

#FFLAGS= -fdefault-real-8 -O -ffree-line-length-256 -frecord-marker=4 $(PIC)



FFLAGS+= -I. -I..

LIBS+=$(FFTW_LIBS)

FILES=$(SRC)/muse_helpers.f90 $(SRC)/mem.f90 $(SRC)/elements.f90 $(SRC)/cooling.f90 \
  $(SRC)/ionize.f90 $(SRC)/molecules.f90 $(SRC)/H2cool.f90 $(SRC)/starprop.f90 \
  $(SRC)/random.f90 $(SRC)/io.f90 $(SRC)/io-old.f $(SRC)/diagnostics.f90 $(SRC)/readparam.f90 \
  $(SRC)/ppread.f90 $(SRC)/blackhole.f90 $(SRC)/stepsystem.f90 $(SRC)/buildnearlist.f90 \
  $(SRC)/pmgrav.f90 $(SRC)/search.f90 $(SRC)/feedback.f90 $(SRC)/fixedhalo.f90 $(SRC)/init.f90 \
  $(SRC)/timestep.f90 $(SRC)/buildtree.f90 $(SRC)/density.f90 $(SRC)/io2.f90 $(SRC)/heating.f90 \
  $(SRC)/entdot.f90 $(SRC)/heco.f90 $(SRC)/ethstep.f90 $(SRC)/stellar.f90 $(SRC)/fcco.f90 \
  $(SRC)/makemapmod.f90 $(SRC)/project.f90 $(SRC)/clean.f90 $(SRC)/sort.f90 $(SRC)/util.f90 $(SRC)/gravity.f90 \
  $(SRC)/tidalsum.f90 $(SRC)/treewalk.f90 $(SRC)/gravsum.f90 $(SRC)/fuvflux.f90

MAPFILES=$(SRC)/util.f90 $(SRC)/makemapmod.f90 $(SRC)/map_helpers.f90

AMUSE_INTERFACE_DEPS=

all: fi_worker

%.o: %.f90 Makefile
	$(MPIFC) $(FFLAGS)  -c -o $@ $< 
	
$(BUILDDIR)/Makefile: 
	mkdir -p $(BUILDDIR)
	ln -s ../src/Makefile  $(BUILDDIR)/Makefile
	
$(BUILDDIR)/libfi.a: $(BUILDDIR)/Makefile $(FILES) Makefile
	mkdir -p $(BUILDDIR)
	$(PYTHON) ./patch_globals_dot_h.py
	make -C $(BUILDDIR) amuse_interface FFLAGS="$(FFLAGS)" FC="$(FC)" VPATH=../src  F77FLAGS="$(F77FLAGS)"

fi_worker: worker_code.f90 $(BUILDDIR)/libfi.a interface.o 
	$(MPIFC) $(FFLAGS) $(SC_FLAGS) $(FS_FLAGS) -I$(BUILDDIR) worker_code.f90 interface.o -o $@ $(BUILDDIR)/libfi.a  $(LIBS) $(SC_FCLIBS) $(FS_LIBS)

fi_worker_gl: worker_code-gl.f90 $(BUILDDIR)/libfi.a interface.o 
	make -C $(BUILDDIR) gl FFLAGS="$(FFLAGS)" FC="$(FC)" GLINCLUDE="$(GLINCLUDE)" GLLIB="$(GLLIB)" VPATH=../src  F77FLAGS="$(F77FLAGS)"
	$(MPIFC) $(FFLAGS) $(SC_FLAGS) $(FS_FLAGS) -I$(BUILDDIR) worker_code-gl.f90 interface.o -o $@  $(GLOBJECTS) \
	$(GLLIB)  $(X11LIB) $(THREADLIB) $(BUILDDIR)/libfi.a $(LIBS) $(SC_FCLIBS) $(FS_LIBS)

$(PERIODIC_BUILDDIR)/Makefile:
	mkdir -p $(PERIODIC_BUILDDIR)
	ln -s ../src/Makefile $(PERIODIC_BUILDDIR)/Makefile
	
$(PERIODIC_BUILDDIR)/libfi.a: $(PERIODIC_BUILDDIR)/Makefile 
	$(PYTHON) ./patch_globals_dot_h.py --periodic
	make -C $(PERIODIC_BUILDDIR) amuse_interface FFLAGS="-I./$(PERIODIC_BUILDDIR) $(FFLAGS)" FC="$(FC)" VPATH=../src   F77FLAGS="$(F77FLAGS)"
	$(PYTHON) ./patch_globals_dot_h.py
	
fi_worker_periodic: worker_code.f90 $(PERIODIC_BUILDDIR)/libfi.a interface.po 
	$(MPIFC) -I./$(PERIODIC_BUILDDIR) $(FFLAGS) $(SC_FLAGS) $(FS_FLAGS) worker_code.f90 interface.po -o $@  $(PERIODIC_BUILDDIR)/libfi.a $(SC_FCLIBS) $(FS_LIBS) $(FFTW_LIBS)

fi_worker_periodic_gl: worker_code-gl.f90 $(PERIODIC_BUILDDIR)/libfi.a interface.po 
	make -C $(PERIODIC_BUILDDIR) gl FFLAGS="$(FFLAGS)" FC="$(FC)" GLINCLUDE="$(GLINCLUDE)" GLLIB="$(GLLIB)" VPATH=../src  F77FLAGS="$(F77FLAGS)"
	$(MPIFC) -I./$(PERIODIC_BUILDDIR) $(FFLAGS) $(SC_FLAGS) $(FS_FLAGS) worker_code-gl.f90 interface.po -o $@  $(PGLOBJECTS) \
	$(GLLIB)  $(X11LIB) $(PTHREADLIB) $(PERIODIC_BUILDDIR)/libfi.a $(LIBS) $(SC_FCLIBS) $(FS_LIBS) $(FFTW_LIBS)

$(OPENMP_BUILDDIR)/Makefile:
	mkdir -p $(OPENMP_BUILDDIR)
	ln -s ../src/Makefile $(OPENMP_BUILDDIR)/Makefile
	
$(OPENMP_BUILDDIR)/libfi.a: $(OPENMP_BUILDDIR)/Makefile 
	make -C $(OPENMP_BUILDDIR) amuse_interface FFLAGS="-I./$(OPENMP_BUILDDIR) $(FFLAGS) $(OPENMP_FCFLAGS)" FC="$(FC)" VPATH=../src   F77FLAGS="$(F77FLAGS) $(OPENMP_FCFLAGS)"

fi_worker_mp: BUILDDIR = $(OPENMP_BUILDDIR)
fi_worker_mp: worker_code.f90 $(OPENMP_BUILDDIR)/libfi.a interface.o 
	$(MPIFC) $(FFLAGS) $(OPENMP_FCFLAGS) $(SC_FLAGS) $(FS_FLAGS) -I$(OPENMP_BUILDDIR) worker_code.f90 interface.o -o $@ $(OPENMP_BUILDDIR)/libfi.a  $(LIBS) $(SC_FCLIBS) $(FS_LIBS)

worker_code.f90: interface.py $(AMUSE_INTERFACE_DEPS)
	$(CODE_GENERATOR) --type=f90 --needs-mpi=false $< FiInterface -o $@

worker_code-gl.f90: interface.py $(AMUSE_INTERFACE_DEPS)
	$(CODE_GENERATOR) --type=f90 $< GlFiInterface -o $@

interface.o:interface.f90
	$(MPIFC) $(FFLAGS) $(SC_FLAGS) -I$(BUILDDIR) -c -o $@ $<
	
interface.po:interface.f90
	$(MPIFC) $(FFLAGS) $(SC_FLAGS) -I./$(PERIODIC_BUILDDIR) -c -o $@ $<


build_map/Makefile: 
	mkdir -p build_map
	ln -s ../src/Makefile  build_map/Makefile

build_map/libmap.a: build_map/Makefile $(MAPFILES) Makefile
	mkdir -p build_map
	make -C build_map libmap.a FFLAGS="$(FFLAGS)" FC="$(FC)" VPATH=../src  F77FLAGS="$(F77FLAGS)"

fi_worker_map: BUILDDIR=./build_map
fi_worker_map: map_worker_code.f90 build_map/libmap.a map_interface.o 
	$(MPIFC) $(FFLAGS) $(SC_FLAGS) $(FS_FLAGS) -I./build_map/ map_worker_code.f90 map_interface.o \
	  -o $@ build_map/libmap.a $(LIBS) $(SC_FCLIBS) $(FS_LIBS)

map_worker_code.f90: interface.py $(AMUSE_INTERFACE_DEPS)
	$(CODE_GENERATOR) --type=f90 --needs-mpi=false $< FiMapInterface -o $@

map_interface.o:map_interface.f90
	$(MPIFC) $(FFLAGS) $(SC_FLAGS) -I$(BUILDDIR) -c -o $@ $<


build_map_mp/Makefile: 
	mkdir -p build_map_mp
	ln -s ../src/Makefile  build_map_mp/Makefile

build_map_mp/libmap.a: build_map_mp/Makefile $(MAPFILES) Makefile
	mkdir -p build_map_mp
	make -C build_map_mp libmap.a FFLAGS="$(FFLAGS) $(OPENMP_FCFLAGS)" FC="$(FC)" VPATH=../src  F77FLAGS="$(F77FLAGS)"

fi_worker_map_mp: BUILDDIR=./build_map_mp
fi_worker_map_mp: map_worker_code.f90 build_map_mp/libmap.a map_interface.o 
	$(MPIFC) $(FFLAGS) $(OPENMP_FCFLAGS) $(SC_FLAGS) $(FS_FLAGS) -I./build_map_mp/ map_worker_code.f90 map_interface.o \
	  -o $@ build_map_mp/libmap.a $(LIBS) $(SC_FCLIBS) $(FS_LIBS)



clean:
	make -C src purge
	$(RM) -f *.o *.po *.pyc *.bck worker worker.f90 glworker glworker.f90
	$(RM) -f worker_code.f90 worker_code_gl.f90 worker_code_sockets.f90 worker_code-gl.f90
	$(RM) -f fi_worker_periodic fi_worker_gl fi_worker fi_worker_mp
	$(RM) -f fi_worker_periodic_gl
	$(RM) -f fi_worker_sockets worker-sockets.f90
	$(RM) -f worker_code_periodic worker_code
	$(RM) -f amuse_worker_module.mod
	$(RM) -f amuseinterface.mod mapinterface.mod
	$(RM) -f forsockets.mod
	$(RM) -f map_worker_code.f90 fi_worker_map fi_worker_map_mp
	$(RM) -Rf $(BUILDDIR)
	$(RM) -Rf $(PERIODIC_BUILDDIR)
	$(RM) -Rf $(OPENMP_BUILDDIR)
	$(RM) -Rf build_map build_map_mp

