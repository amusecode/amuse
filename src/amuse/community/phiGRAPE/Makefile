# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(strip $(AMUSE_DIR)),)
    AMUSE_DIR = ../../../..
endif
-include ${AMUSE_DIR}/config.mk


ifeq ($(strip $(MPICXX)),)
    MPICXX = mpicxx
endif
ifeq ($(strip $(MPICC)),)
    MPICC = mpicc
endif
ifeq ($(strip $(MPIFC)),)
    MPIFC = mpif90
endif
ifeq ($(strip $(MPIF90)),)
    MPIF90 = $(MPIFC)
endif
FC = $(MPIF90)
CC = $(MPICC)
ifeq ($(strip $(GL_PATH)),)
    GL_PATH = /data/pelupes/libraries/f90gl-1.2.11
endif


CODE_GENERATOR = $(AMUSE_DIR)/build.py

SRCDIR = src

MPI_BUILDDIR=build_mpi
BUILDDIR = build

CODELIB=$(BUILDDIR)/libphigrape.a

GLOBJECTS=src/fthread.o src/viewer.o src/start_viewer.o src/libpt.a

ifeq ($(strip $(CUDA_SDK)),)
    CUDA_SDK = /disks/koppoel1/CUDA23/cudasdk/C
endif
ifeq ($(strip $(CUDA_TK)),)
    CUDA_TK = /disks/koppoel1/CUDA23/cuda
endif

ifeq ($(strip $(CUDA_LIBDIRS)),)
    CUDA_LIBDIRS = -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64
endif
ifeq ($(strip $(CUDA_LIBS)),)
    CUDA_LIBS = -lcudart
endif

ifeq ($(strip $(SAPPORO_LIBDIRS)),)
    SAPPORO_LIBDIRS = -L$(AMUSE_DIR)/lib/sapporo_light
endif
ifeq ($(strip $(OPENMP_CFLAGS)),)
    OPENMP_CFLAGS = SAPPORO_LIBS ?= $(SAPPORO_LIBDIRS) -lsapporo $(OPENMP_CFLAGS)
endif
ifeq ($(strip $(GRAPE6_LIBDIRS)),)
    GRAPE6_LIBDIRS = -L$(GRAPE6)/lib
endif

BOOSTLIBS = -L/opt/local/lib -lboost_thread-mt
LLIBSDIR =  $(GRAPE6_LIBDIRS)
GPULIBS  =  $(SAPPORO_LIBS)  $(CUDA_LIBDIRS) $(CUDA_LIBS) -lpthread -L/usr/lib64 -lstdc++
ifeq ($(strip $(GRAPE6LIBS)),)
    GRAPE6LIBS = -lg6lx
endif
LIBS     +=  -lm -lpthread

LIBRARIES= $(LLIBSDIR) $(LIBS)

GLINCLUDE = -I/usr/include/GL -I$(GL_PATH)/include/GL/
GLLIB = -L$(GL_PATH)/lib -lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut

ifeq ($(strip $(X11LIB)),)
    X11LIB = -L/usr/X11R6/lib64 -lXaw -lXt -lXmu -lXi -lXext -lX11
endif

ifeq ($(strip $(THREADLIB)),)
    THREADLIB = -L. -lpthread src/libpt.a
endif

ifeq ($(strip $(PGLIBS)),)
    PGLIBS = -L. -lpg5
endif

ifeq ($(strip $(G6LIBS)),)
    G6LIBS = -L$(AMUSE_DIR)/lib/g6 -lg6
endif

CFLAGS += -g  $(OPTS)
FFLAGS += -g $(OPTS)
FFLAGS += -I./src/

ifeq ($(findstring gfortran, $(shell $(FC) --help)), gfortran)
else
GLLIB += -lguide
endif

all: debuginfo phigrape_worker 

debuginfo:
	@echo "making phigrape..; some debugging info follows:"
	@echo $(SC_MPI_CLIBS)


$(BUILDDIR)/Makefile: $(SRCDIR)/Makefile 
	-mkdir $(BUILDDIR) 
	cp $(SRCDIR)/Makefile $(BUILDDIR)/Makefile
	
$(CODELIB): $(BUILDDIR)/Makefile
	make -C $(BUILDDIR) all "FC=$(FC)" "OPTS=-DNOMPI" "CC=$(CC)" VPATH=../src

phigrape_worker: worker_code.f90 interface.o $(CODELIB)
	$(MPIF90) $(FFLAGS) $(SC_FLAGS) $(FS_FLAGS) $(LDFLAGS) -DNOMPI $^ $(CODELIB) -o $@ $(OBJECTS) $(G6LIBS) $(LIBRARIES) $(SC_MPI_CLIBS) $(FS_LIBS)  $(LIBS)

$(MPI_BUILDDIR)/Makefile: $(SRCDIR)/Makefile 
	-mkdir $(MPI_BUILDDIR) 
	cp $(SRCDIR)/Makefile $(MPI_BUILDDIR)/Makefile
	
$(MPI_BUILDDIR)/libphigrape.a:  $(MPI_BUILDDIR)/Makefile
	make -C $(MPI_BUILDDIR) all "FC=$(FC)" "OPTS=-DMPI" "CC=$(CC)" VPATH=../src
	
phigrape_worker_mpi: worker_code.f90 interface.mpio $(MPI_BUILDDIR)/libphigrape.a
	$(MPIF90) $(FFLAGS) $(SC_FLAGS) $(FS_FLAGS) $(LDFLAGS) -DMPI $^ $(MPI_BUILDDIR)/libphigrape.a -o $@ $(OBJECTS) $(G6LIBS) $(LIBRARIES) $(SC_MPI_CLIBS) $(FS_LIBS)  $(LIBS)

phigrape_worker_gl: glworker_code.f90 interface.o $(CODELIB)
	make -C src gl "FC=$(FC)" "OPTS=-DNOMPI" "CC=$(CC)"
	$(MPIF90) $(FFLAGS) -DNOMPI $(SC_FLAGS) $(FS_FLAGS) $(LDFLAGS) $^ $(CODELIB) -o $@ $(OBJECTS) $(G6LIBS) $(LIBRARIES) \
	$(GLOBJECTS) $(GLLIB)  $(X11LIB) $(THREADLIB) $(LIBS) $(SC_MPI_CLIBS) $(FS_LIBS)  $(LIBS)

phigrape_worker_gpu: worker_code.f90 interface.o $(CODELIB)
	$(MPIF90) $(FFLAGS) -DNOMPI $(SC_FLAGS) $(FS_FLAGS) $(LDFLAGS) $^ $(CODELIB) -o $@ $(OBJECTS) $(LIBRARIES) $(GPULIBS) $(SC_MPI_CLIBS) $(FS_LIBS)  $(LIBS)

phigrape_worker_gl_gpu: glworker_code.f90 interface.o $(CODELIB)
	make -C src gl "FC=$(FC)" "OPTS=-DNOMPI" "CC=$(CC)"
	$(MPIF90) $(FFLAGS) -DNOMPI $(SC_FLAGS) $(FS_FLAGS) $(LDFLAGS) $^ $(CODELIB) -o $@ $(OBJECTS) $(LIBRARIES) \
	$(GLOBJECTS) $(GLLIB)  $(X11LIB) $(THREADLIB) $(LIBS) $(GPULIBS) $(SC_MPI_CLIBS) $(FS_LIBS)  $(LIBS)

phigrape_worker_grape: worker_code.f90 interface.o $(CODELIB)
	$(MPIF90) $(FFLAGS) -DNOMPI $(SC_FLAGS) $(FS_FLAGS) $(LDFLAGS) $^ $(CODELIB) -o $@ $(OBJECTS) $(LIBRARIES) $(GRAPE6LIBS) $(SC_FCLIBS) $(FS_LIBS)  $(LIBS)

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py PhiGRAPEInterface -o $@

glworker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py PhiGRAPEInterfaceGL -o $@

clean:
	make -C src clean
	rm -f *.pyc *.o *.mpio
	rm -Rf build
	rm -Rf build_mpi
	rm -f worker_code_gpu worker_code_grape
	rm -f worker_code worker_code.f90
	rm -f worker_code-sockets.f90 phigrape_worker_sockets
	rm -f glworker_code glworker_code.f90
	rm -f glworker_code_gpu glworker_code_grape
	rm -f phigrape_worker phigrape_worker_gl phigrape_worker_gl_gpu phigrape_worker_gpu
	rm -f phigrape_worker_mpi 
	rm -f phigrape_worker_gpu_sockets
	rm -f forsockets.mod

__init__.py:
	touch $@

%.o: %.f90
	$(FC)  $(FFLAGS) -c $<

%.mpio: %.F
	$(FC)  $(SC_FLAGS) -DMPI $(FFLAGS) -c $< -o $@

%.o: %.F
	$(FC)  $(SC_FLAGS) -DNOMPI $(FFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) $(SC_FLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

.F.o:
	$(FC) $(FFLAGS) $(SC_FLAGS) -c $<
