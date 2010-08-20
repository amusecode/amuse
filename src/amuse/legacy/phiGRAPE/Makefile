FC = mpif90
CC = mpicc

OPTS+=-DNOMPI

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

SRCDIR = src

CODELIB=$(SRCDIR)/libphigrape.a

GLOBJECTS=src/fthread.o src/viewer.o src/start_viewer.o src/libpt.a

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

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

GLINCLUDE = -I/usr/include/GL \
-I/data/pelupes/libraries/f90gl-1.2.11/include/GL/

GLLIB = -L/data/pelupes/libraries/f90gl-1.2.11/lib \
-lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut  -lguide 

X11LIB ?= -L/usr/X11R6/lib64 -lXaw -lXt -lXmu -lXi -lXext -lX11

THREADLIB ?= -L. -lpthread src/libpt.a

PGLIBS ?= -L. -lpg5

G6LIBS ?= -L$(AMUSE_DIR)/lib/g6 -lg6

CFLAGS += -g  $(OPTS)
FFLAGS += -g $(OPTS)
FFLAGS += -I./src/


all: debuginfo worker_code 

debuginfo:
	@echo "making phigrape..; some debugging info follows:"
	@echo $(SCLIBS)

$(CODELIB):
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"

worker_code: worker_code.f90 phigrape_code.o $(CODELIB)
	mpif90 $(FFLAGS) $(SCINC) $^ $(CODELIB) -o $@ $(OBJECTS) $(G6LIBS) $(LIBRARIES) $(SCLIBS)

glworker_code: glworker_code.f90 phigrape_code.o $(CODELIB)
	make -C src gl "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $(FFLAGS) $^ $(CODELIB) -o $@ $(OBJECTS) $(G6LIBS) $(LIBRARIES) \
	$(GLOBJECTS) $(GLLIB)  $(X11LIB) $(THREADLIB) $(LIBS)

worker_code_gpu: worker_code.f90 phigrape_code.o $(CODELIB)
	mpif90 $(FFLAGS) $^ $(CODELIB) -o $@ $(OBJECTS) $(LIBRARIES) $(GPULIBS)

glworker_code_gpu: glworker_code.f90 phigrape_code.o $(CODELIB)
	make -C src gl "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $(FFLAGS) $^ $(CODELIB) -o $@ $(OBJECTS) $(LIBRARIES) \
	$(GLOBJECTS) $(GLLIB)  $(X11LIB) $(THREADLIB) $(LIBS) $(GPULIBS)


worker_code_grape: worker_code.f90 phigrape_code.o $(CODELIB)
	mpif90 $(FFLAGS) $^ $(CODELIB) -o $@ $(OBJECTS) $(LIBRARIES) $(GRAPE6LIBS)

worker_code_phantom_grape: worker_code.f90 phigrape_code.o $(CODELIB)
	mpif90 $(FFLAGS) $^ $(CODELIB) -o $@ $(OBJECTS) $(LIBRARIES) $(PGLIBS)

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py PhiGRAPEInterface -o $@

glworker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py PhiGRAPEInterfaceGL -o $@

clean:
	make -C src clean
	rm -f *.pyc *.o
	rm -f worker_code_gpu worker_code_grape
	rm -f worker_code worker_code.f90
	rm -f glworker_code glworker_code.f90
	rm -f glworker_code_gpu glworker_code_grape

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
