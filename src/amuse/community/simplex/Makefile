MPICXX ?= mpicxx
export MPICXX

AMUSE_DIR?=../../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py

SRCDIR = src
OBJDIR = $(SRCDIR)/obj
APP		 = $(SRCDIR)/bin/SimpleX
PLUGINS = $(SRCDIR)/plugins

OBJS = $(OBJDIR)/amuse_interface.o $(OBJDIR)/SimpleX.o \
       $(OBJDIR)/Structs.o $(OBJDIR)/Common.o $(OBJDIR)/h5w_serial.o \
       $(OBJDIR)/keyValue.o $(OBJDIR)/tree_structures.o $(OBJDIR)/hilbert.o 
GLOBJS = $(OBJDIR)/amuse_interface_gl.o

HEALPIX_PATH = $(PLUGINS)/HEALPix/
QHULL_PATH = $(PLUGINS)/QHull/

HDF5_PATH = /Users/jppaarde/sw/hdf5-1.6.10/installed
GSL_PATH = /usr/local/scisoft

LIBS =  -L${HDF5_PATH}/lib -lhdf5 -lz \
        -L${QHULL_PATH} -lqhull \
        -L${HEALPIX_PATH}/lib  -lhealpix_cxx \
        -L${GSL_PATH}/lib -lgslcblas -lgsl 

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

worker.cc: interface.py
	$(CODE_GENERATOR) --type=c $< SimpleXInterface -o $@

worker.h: interface.py
	$(CODE_GENERATOR) --type=H $< SimpleXInterface -o $@

worker_gl.cc: interface.py
	$(CODE_GENERATOR) --type=c $< GlSimpleXInterface -o $@

worker_gl.h: interface.py
	$(CODE_GENERATOR) --type=H $< GlSimpleXInterface -o $@

simplex_worker: worker.cc worker.h $(OBJS)
	$(MPICXX) $^ -o $@ $(LIBS)

simplex_worker_gl: worker_gl.cc worker_gl.h $(OBJS) $(GLOBJS)
	$(MPICXX) $^ -o $@ $(LIBS) $(GLLIBS)

$(OBJDIR)/amuse_interface.o: $(SRCDIR)/src/amuse_interface.cpp $(SRCDIR)/src/amuse_interface.h
	cd $(SRCDIR) && $(MAKE) amuse_interface
$(OBJDIR)/amuse_interface_gl.o: $(SRCDIR)/src/amuse_interface_gl.cpp
	cd $(SRCDIR) && $(MAKE) amuse_interface_gl
