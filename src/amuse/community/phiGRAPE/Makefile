MPICXX ?= mpicxx
MPICC ?= mpicc
MPIF90 ?= mpif90
FC = $(MPIF90)
CC = $(MPICC)
OPTS+=-DNOMPI
GL_PATH ?= /data/pelupes/libraries/f90gl-1.2.11

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

SRCDIR = src

CODELIB=$(SRCDIR)/libphigrape.a

GLOBJECTS=src/fthread.o src/viewer.o src/start_viewer.o src/libpt.a

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

FSLIBS ?= -L$(AMUSE_DIR)/lib/forsockets -lforsockets
FSINC ?= -I$(AMUSE_DIR)/lib/forsockets

CUDA_SDK ?= /disks/koppoel1/CUDA23/cudasdk/C
CUDA_TK  ?= /disks/koppoel1/CUDA23/cuda

CUDA_LIBDIRS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64
CUDA_LIBS ?= -lcudart

SAPPORO_LIBDIRS ?= -L$(AMUSE_DIR)/lib/sapporo_light
GRAPE6_LIBDIRS ?= -L$(GRAPE6)/lib

BOOSTLIBS = -L/opt/local/lib -lboost_thread-mt
LLIBSDIR = $(SAPPORO_LIBDIRS) $(CUDA_LIBDIRS) $(GRAPE6_LIBDIRS)
GPULIBS  =  $(CUDA_LIBS) -lsapporo  -lpthread -L/usr/lib64 -lstdc++
GRAPE6LIBS  ?= -lg6lx
LIBS     =  -lm -lpthread

LIBRARIES= $(LLIBSDIR) $(LIBS)

GLINCLUDE = -I/usr/include/GL -I$(GL_PATH)/include/GL/
GLLIB = -L$(GL_PATH)/lib -lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut

X11LIB ?= -L/usr/X11R6/lib64 -lXaw -lXt -lXmu -lXi -lXext -lX11

THREADLIB ?= -L. -lpthread src/libpt.a

PGLIBS ?= -L. -lpg5

G6LIBS ?= -L$(AMUSE_DIR)/lib/g6 -lg6

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
	@echo $(SCLIBS)

$(CODELIB):
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"

phigrape_worker: worker_code.f90 phigrape_code.o $(CODELIB)
	$(MPIF90) $(FFLAGS) $(SCINC) $^ $(CODELIB) -o $@ $(OBJECTS) $(G6LIBS) $(LIBRARIES) $(SCLIBS)

phigrape_worker_sockets: worker_code-sockets.f90 phigrape_code.o $(CODELIB)
	$(MPIF90) $(FFLAGS) $(SCINC) $(FSINC) $^ $(CODELIB) -o $@ $(OBJECTS) $(G6LIBS) $(LIBRARIES) $(SCLIBS) $(FSLIBS)

phigrape_worker_gl: glworker_code.f90 phigrape_code.o $(CODELIB)
	make -C src gl "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	$(MPIF90) $(FFLAGS) $(SCINC) $^ $(CODELIB) -o $@ $(OBJECTS) $(G6LIBS) $(LIBRARIES) \
	$(GLOBJECTS) $(GLLIB)  $(X11LIB) $(THREADLIB) $(LIBS) $(SCLIBS)

phigrape_worker_gpu: worker_code.f90 phigrape_code.o $(CODELIB)
	$(MPIF90) $(FFLAGS) $(SCINC) $^ $(CODELIB) -o $@ $(OBJECTS) $(LIBRARIES) $(GPULIBS) $(SCLIBS)

phigrape_worker_gpu_sockets: worker_code-sockets.f90 phigrape_code.o $(CODELIB)
	$(MPIF90) $(FFLAGS) $(SCINC) $(FSINC) $^ $(CODELIB) -o $@ $(OBJECTS) $(LIBRARIES) $(GPULIBS) $(SCLIBS) $(FSLIBS)

phigrape_worker_gl_gpu: glworker_code.f90 phigrape_code.o $(CODELIB)
	make -C src gl "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	$(MPIF90) $(FFLAGS) $(SCINC) $^ $(CODELIB) -o $@ $(OBJECTS) $(LIBRARIES) \
	$(GLOBJECTS) $(GLLIB)  $(X11LIB) $(THREADLIB) $(LIBS) $(GPULIBS) $(SCLIBS)

phigrape_worker_grape: worker_code.f90 phigrape_code.o $(CODELIB)
	$(MPIF90) $(FFLAGS) $^ $(CODELIB) -o $@ $(OBJECTS) $(LIBRARIES) $(GRAPE6LIBS) $(SCLIBS)

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py PhiGRAPEInterface -o $@

worker_code-sockets.f90: interface.py
	$(CODE_GENERATOR) --type=f90 --mode=sockets interface.py PhiGRAPEInterface -o $@

glworker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py PhiGRAPEInterfaceGL -o $@

clean:
	make -C src clean
	rm -f *.pyc *.o
	rm -f worker_code_gpu worker_code_grape
	rm -f worker_code worker_code.f90
	rm -f worker_code-sockets.f90 phigrape_worker_sockets
	rm -f glworker_code glworker_code.f90
	rm -f glworker_code_gpu glworker_code_grape
	rm -f phigrape_worker phigrape_worker_gl phigrape_worker_gl_gpu phigrape_worker_gpu

__init__.py:
	touch $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) $(SCINC) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

.F.o:
	$(FC) $(FFLAGS) $(SCINC) -c $<
