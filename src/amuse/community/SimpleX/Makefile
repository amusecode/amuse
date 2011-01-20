SRCDIR = src
OBJDIR = $(SRCDIR)/obj
APP		 = $(SRCDIR)/bin/SimpleX

OBJS = $(OBJDIR)/amuse_interface.o $(OBJDIR)/SimpleX.o \
       $(OBJDIR)/Structs.o $(OBJDIR)/Common.o $(OBJDIR)/h5w_serial.o \
       $(OBJDIR)/keyValue.o $(OBJDIR)/tree_structures.o $(OBJDIR)/hilbert.o 
GLOBJS = $(OBJDIR)/amuse_interface_gl.o

PLUGINS = /data2/pelupes/amuse/trunk/src/amuse/community/SimpleX/src/plugins
HDF5 = /data1/pelupes/libraries/hdf5-1.6.9/
HEALPIX_PATH = $(PLUGINS)/HEALPix/Healpix_2.11/src/cxx/generic_gcc/
GSL_PATH = /usr/local/scisoft

LIBS =  -L${HDF5}/lib -lhdf5 -lz \
        -L${PLUGINS}/QHull/src -lqhull \
        -L${HEALPIX_PATH}/lib  -lhealpix_cxx \
        -L${GSL_PATH}/lib -lgslcblas -lgsl 

GLLIBS=-lGL -lGLU -L/home/inti/libraries/freeglut-2.4.0/src/.libs -lglut \
-L/usr/X11R6/lib -lXaw -lXt -lXmu -lXi -lXext -lX11

#AMUSE_DIR?=../../../..
AMUSE_DIR?=/Users/jppaarde/work/code/amuse-svn

CODE_GENERATOR = $(AMUSE_DIR)/build.py

#all:	compile muse_worker muse_worker_gl
all:	compile muse_worker

compile: $(APP)

$(APP):
	cd $(SRCDIR) && $(MAKE)
	
clean:
	rm -f *.so *.o *.pyc muse_worker.cc muse_worker.h muse_worker_gl.cc muse_worker_gl.h
	$(RM) muse_worker muse_worker_gl *~
	cd $(SRCDIR) && $(MAKE) clean

muse_worker.cc: amuse_simplex.py
	$(CODE_GENERATOR) --type=c amuse_simplex.py SimpleX -o $@

muse_worker.h: amuse_simplex.py
	$(CODE_GENERATOR) --type=h amuse_simplex.py SimpleX -o $@

muse_worker_gl.cc: amuse_simplex_gl.py
	$(CODE_GENERATOR) --type=c amuse_simplex_gl.py SimpleX -o $@

muse_worker_gl.h: amuse_simplex_gl.py
	$(CODE_GENERATOR) --type=h amuse_simplex_gl.py SimpleX -o $@

muse_worker: muse_worker.cc muse_worker.h $(OBJS) $(SRCDIR)/src/amuse_interface.cpp 
	cd $(SRCDIR) && $(MAKE) amuse_interface
	mpicxx $@.cc $(OBJS) -o $@ $(LIBS) $(HEALPIX_PATH)/lib/libhealpix_cxx.a 

muse_worker_gl: muse_worker_gl.cc muse_worker_gl.h $(OBJS) $(SRCDIR)/src/amuse_interface.cpp $(SRCDIR)/src/amuse_interface_gl.cpp
	cd $(SRCDIR) && $(MAKE) amuse_interface
	cd $(SRCDIR) && $(MAKE) amuse_interface_gl
	mpicxx $(GLOBJS) $@.cc $(LIBS) $(GLLIBS) $(OBJS)  $(HEALPIX_PATH)/lib/libhealpix_cxx.a -o $@

$(SRCDIR)/obj/amuse_interface.o:
	cd $(SRCDIR) && $(MAKE) amuse_interface
$(SRCDIR)/obj/amuse_interface_gl.o:
	cd $(SRCDIR) && $(MAKE) amuse_interface_gl


.cc.o: $<
	g++ $(CXXFLAGS) -c -o $@ $< 
.C.o: $<
	g++ $(CXXFLAGS) -c -o $@ $<
