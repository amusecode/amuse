MPICXX ?= mpicxx

AMUSE_DIR?=../../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py

SRCDIR = src
OBJDIR = $(SRCDIR)/obj
APP		 = $(SRCDIR)/bin/SimpleX
PLUGINS = $(SRCDIR)/plugins

HDF5_INCL = -DH5_USE_16_API
HDF5_LIB =  -lhdf5 -lz 

OBJS = $(OBJDIR)/amuse_interface.o $(OBJDIR)/SimpleX.o \
       $(OBJDIR)/Structs.o $(OBJDIR)/Common.o $(OBJDIR)/h5w_serial.o \
       $(OBJDIR)/keyValue.o $(OBJDIR)/tree_structures.o $(OBJDIR)/hilbert.o 
GLOBJS = $(OBJDIR)/amuse_interface_gl.o

HEALPIX_PATH = $(PLUGINS)/HEALPix/Healpix_2.11/src/cxx/generic_gcc/

LIBS =  ${HDF5_LIB} \
        -L${PLUGINS}/QHull/src -lqhull \
        -L${HEALPIX_PATH}/lib  -lhealpix_cxx \
        -lgslcblas -lgsl 

GL_PATH ?= /home/inti/libraries/f90gl-1.2.11
GLINCLUDE = -I/usr/include/GL -I$(GL_PATH)/include/GL/
GLLIB = -L$(GL_PATH)/lib -lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut  #-lguide 
GLLIBS= $(GLLIB) -lGL -lGLU -L/home/inti/libraries/freeglut-2.4.0/src/.libs -lglut \
-L/usr/X11R6/lib -lXaw -lXt -lXmu -lXi -lXext -lX11 -ltiff

all:	compile simplex_worker

compile: $(APP)

$(APP):
	cd $(SRCDIR) && $(MAKE)
	
clean:
	$(RM) *.so *.o *.pyc worker.cc worker.h worker_gl.cc worker_gl.h
	$(RM) simplex_worker simplex_worker_gl *~
	cd $(SRCDIR) && $(MAKE) clean

worker.cc: amuse_simplex.py
	$(CODE_GENERATOR) --type=c amuse_simplex.py SimpleX -o $@

worker.h: amuse_simplex.py
	$(CODE_GENERATOR) --type=h amuse_simplex.py SimpleX -o $@

worker_gl.cc: amuse_simplex_gl.py
	$(CODE_GENERATOR) --type=c amuse_simplex_gl.py SimpleX -o $@

worker_gl.h: amuse_simplex_gl.py
	$(CODE_GENERATOR) --type=h amuse_simplex_gl.py SimpleX -o $@

simplex_worker: worker.cc worker.h $(OBJS)
	$(MPICXX) $^ -o $@ $(LIBS)

simplex_worker_gl: worker_gl.cc worker_gl.h $(OBJS) $(GLOBJS)
	$(MPICXX) $^ -o $@ $(LIBS) $(GLLIBS)

$(OBJDIR)/amuse_interface.o: $(SRCDIR)/src/amuse_interface.cpp
	cd $(SRCDIR) && $(MAKE) amuse_interface
$(OBJDIR)/amuse_interface_gl.o: $(SRCDIR)/src/amuse_interface_gl.cpp
	cd $(SRCDIR) && $(MAKE) amuse_interface_gl
