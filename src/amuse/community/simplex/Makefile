# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

MPICXX ?= mpicxx
export MPICXX


SRCDIR = src
OBJDIR = $(SRCDIR)/obj
APP		 = $(SRCDIR)/bin/SimpleX
PLUGINS = $(SRCDIR)/plugins

OBJS = $(OBJDIR)/amuse_interface.o $(OBJDIR)/SimpleX.o \
       $(OBJDIR)/Structs.o $(OBJDIR)/rates.o $(OBJDIR)/Common.o $(OBJDIR)/h5w_serial.o \
       $(OBJDIR)/keyValue.o $(OBJDIR)/tree_structures.o $(OBJDIR)/hilbert.o 

GLOBJS = $(OBJDIR)/amuse_interface_gl.o

HEALPIX_PATH = $(PLUGINS)/HEALPix/
QHULL_PATH = $(PLUGINS)/qhull-2012.1/

HDF5_FLAGS ?= -I$(HDF5_PATH)/include
HDF5_LIBS  ?= -L$(HDF5_PATH)/lib -lhdf5 -lz

GSL_PATH = /usr/local/scisoft

# here we use ?= to allow amuse configure to override these
# variables
GSL_FLAGS ?= -I$(GSL_PATH)/include
GSL_LIBS ?= -L$(GSL_PATH)/lib -lgsl -lgslcblas

LIBS =  $(HDF5_LIBS) \
        $(GSL_LIBS) \
        -L$(QHULL_PATH)/lib -lqhullstatic \
        -L$(HEALPIX_PATH)/lib  -lhealpix_cxx

GLINCLUDE = -I/usr/include/GL -I$(GL_PATH)/include/GL/
GLLIB = -L$(GL_PATH)/lib -lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut  #-lguide 
GLLIBS= $(GLLIB) -lGL -lGLU -L/home/inti/libraries/freeglut-2.4.0/src/.libs -lglut \
-L/usr/X11R6/lib -lXaw -lXt -lXmu -lXi -lXext -lX11 -ltiff

all:	compile simplex_worker

compile: $(APP)

$(APP):
	cd $(SRCDIR) && $(MAKE) CC="$(CC) $(CFLAGS)" CXX="$(CXX) $(CXXFLAGS)"

clean:
	$(RM) *.so *.o *.pyc worker.cc worker.h worker_gl.cc worker_gl.h
	$(RM) worker_code.cc worker_code.h worker_code_gl.cc worker_code_gl.h
	$(RM) simplex_worker simplex_worker_gl *~
	cd $(SRCDIR) && $(MAKE) clean

distclean: clean
	cd $(SRCDIR) && $(MAKE) veryclean
	
worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c $< SimpleXInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h $< SimpleXInterface -o $@

worker_code-gl.cc: interface.py
	$(CODE_GENERATOR) --type=c $< GlSimpleXInterface -o $@

worker_code-gl.h: interface.py
	$(CODE_GENERATOR) --type=h $< GlSimpleXInterface -o $@

simplex_worker: worker_code.cc worker_code.h $(OBJS)
	$(MPICXX) $(LDFLAGS) $< $(OBJS) -o $@ $(LIBS)

#simplex_worker_gl: worker_code-gl.cc worker_code-gl.h $(OBJS) $(GLOBJS)
#	$(MPICXX) $^ -o $@ $(LIBS) $(GLLIBS)

$(OBJDIR)/amuse_interface.o: $(SRCDIR)/src/amuse_interface.cpp $(SRCDIR)/src/amuse_interface.h
	cd $(SRCDIR) && $(MAKE) amuse_interface
	
$(OBJDIR)/amuse_interface_gl.o: $(SRCDIR)/src/amuse_interface_gl.cpp
	cd $(SRCDIR) && $(MAKE) amuse_interface_gl
