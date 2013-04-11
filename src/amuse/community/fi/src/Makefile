## makefile for treesph code
AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcond
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

FC=ifort
F2PY=f2py
PIC=

OPENMP_FCFLAGS ?=  -openmp -openmp_report0 

# gfortran flags
# FFLAGS2=-fdefault-real-8 -ffree-line-length-256  -frecord-marker=4

# ifort flags

FFLAGS2= -r8 -O2 -ip -u  $(PIC)

FFLAGS= $(FFLAGS2) $(OPENMP_FCFLAGS)

APPS= fi glfi snapviewer hicube map sfh sfh2 rgabe togabe \
      extract rsimple rgadget combine_simple snapreader \
      tosimplex tosimple

# -----------------------------------------------------------------------
COMMONFILES = mem.o elements.o cooling.o ionize.o molecules.o H2cool.o starprop.o \
  random.o io-old.o diagnostics.o readparam.o ppread.o \
  blackhole.o stepsystem.o buildnearlist.o pmgrav.o \
  search.o feedback.o fixedhalo.o init.o timestep.o buildtree.o \
  density.o heating.o entdot.o heco.o ethstep.o \
  stellar.o fcco.o makemapmod.o project.o clean.o sort.o util.o \
  gravity.o tidalsum.o treewalk.o gravsum.o fuvflux.o

LIBFILES = $(COMMONFILES) muse_io2_stub.o muse_io_stub.o hash.o muse_helpers.o

MAPLIBFILES = makemapmod.o hash.o map_helpers.o util.o

FILES = $(COMMONFILES) io.o io2.o

SFILES=extern_dummy.o 

GL_PATH ?= /home/inti/libraries/f90gl-1.2.11
GLINCLUDE = -I/usr/include/GL -I$(GL_PATH)/include/GL/
GLLIB = -L$(GL_PATH)/lib -lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut 

X11LIB = # -L/usr/lib64 -lXaw -lXt -lXmu -lXi -lXext -lX11

FITSLIB= -lcfitsio

FFTWLIB =  -lfftw3 -lfftw3_threads

LIBS= $(FFTWLIB) $(FITSLIB)

THREAD_OPTIONS = -L. -lpthread -lpt

ELIBS=-lguide /opt/intel/fc/9.0/lib/libsvml.a

F2PYFLAGS=--quiet --f90exec=$(FC) --f77exec=$(FC) \
 --f90flags=\"$(FFLAGS)\" --f77flags=\"$(FFLAGS)\"

# ------------------------------------------------------------------------------

all: $(APPS)

%.o: %.f90 globals.h Makefile
	$(FC) $(FFLAGS) $(SCINC)  -c -o $@ $< 

%.o: %.f globals.h Makefile
	$(FC) $(F77FLAGS) $(FFLAGS) $(SCINC) -c -o $@ $< 

stars3d.o: stars3d.f90 mem.o
	$(FC) $(FFLAGS) $(GLINCLUDE) -c -o $@ $<

fi_viewer.o: fi_viewer.f90 mem.o
	$(FC) $(FFLAGS) $(GLINCLUDE) -c -o $@ $<

snapviewer.o: snapviewer.f90
	$(FC) $(FFLAGS) $(GLINCLUDE) -c -o $@ $<

fi: $(FILES) $(SFILES) fi.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

glfi: $(FILES) $(SFILES) glfi.o fthread.o stars3d.o libpt.a
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) $(GLLIB) $(X11_OPTIONS) $(THREAD_OPTIONS) $(X11LIB)

snapviewer: $(FILES) $(SFILES) snapviewer.o 
	$(FC) $(FFLAGS) -o $@ $^ $(GLLIB) $(X11LIB) $(LIBS)

hicube: $(FILES) $(SFILES) hicube.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

map: $(FILES) $(SFILES) makemap.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

extract: $(FILES) $(SFILES) extract.o varia.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

sfh: sfh.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

sfh2: sfh2.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

rgadget: mem.o io2.o  rgadget.o 
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

rsimple: mem.o io2.o rsimple.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

tosimple: mem.o io2.o tosimple.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

combine_simple: combine_simple.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

snapreader: snapreader.o $(FILES) $(SFILES)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

togabe: $(FILES) $(SFILES) togabe.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

togabe2: $(FILES) $(SFILES) sphray_io.o togabe2.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

tosimplex: $(FILES) $(SFILES) tosimplex.o varia.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

rgabe: rgabe.o io2.o mem.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)
rgabe2: sphray_io.o rgabe2.o io2.o mem.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

muse: muse_dynamics_interface muse_dynamics_interface_gl \
      muse_sph_interface muse_sph_interface_gl
	rm -f *.so

amuse_interface: libfi.a

libfi.a: $(LIBFILES)
	ar crs $@ $^

libmap.a: $(MAPLIBFILES)
	ar crs $@ $^

muse_dynamics_interface: muse_interface.f90 $(FILES) hash.o muse_helpers.o
	@echo $(F2PY) -c $(F2PYFLAGS) \
         -m $@#  $^ $(LIBS) $(ELIBS) > muse_dynamics_interface.txt

muse_dynamics_interface_gl: muse_interface.f90 muse_interface_gl.f90 \
                            $(FILES) hash.o muse_helpers.o fthread.o stars3d.o libpt.a
	@echo $(F2PY) -c $(F2PYFLAGS) -m $@#  $^ \
            $(LIBS) $(ELIBS) $(GLLIB) $(X11_OPTIONS) \
            $(THREAD_OPTIONS) $(X11LIB) > muse_dynamics_interface_gl.txt

muse_sph_interface: muse_interface.f90 $(FILES) hash.o muse_helpers.o
#	$(F2PY) -c $(F2PYFLAGS) -m $@0  $^ $(LIBS) $(ELIBS)
	@echo $(F2PY) -c $(F2PYFLAGS) -m $@#  $^ \
            $(LIBS) $(ELIBS) > muse_sph_interface.txt

muse_sph_interface_gl: muse_interface.f90 muse_interface_gl.f90 \
                       $(FILES) hash.o muse_helpers.o fthread.o stars3d.o libpt.a
	@echo $(F2PY) -c $(F2PYFLAGS) -m $@#  $^ \
            $(LIBS) $(ELIBS) $(GLLIB) $(X11_OPTIONS) \
            $(THREAD_OPTIONS) $(X11LIB) > muse_sph_interface_gl.txt

objects: $(FILES) $(SFILES)

libpt.a: pt.o ptf77.o pt.h
	ar crs libpt.a pt.o ptf77.o

gl: fi_viewer.o libpt.a fthread.o

clean: 
	rm -f *.o *.mod libpt.a libfi.a *.bck *.so *.pyc 
	rm -f muse_dynamics_interface.txt muse_dynamics_interface_gl.txt
	rm -f muse_sph_interface.txt muse_sph_interface_gl.txt

purge: clean
	rm -f $(APPS)

install: $(APPS) clean
	cp $(APPS) ../bin
	rm -f $(APPS)

