
FC = mpif90
CC = mpicc

OPTS+=-DNOMPI

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

SRCDIR = src

OBJECTS = \
 src/allgather_ap.o src/cmcorr.o src/corrector.o src/cputime.o \
 src/energy.o src/get_min_t.o src/gravity.o src/initgrape.o \
 src/initmpi.o src/initpars.o src/output.o \
 src/predictor.o src/readbodydata.o src/readpars.o src/readrestartfile.o \
 src/selectactive.o src/selectsync.o src/sendbodies2grape.o \
 src/sendbodydata.o \
 src/sendinpdata.o src/sumforces.o src/sync.o \
 src/timestep.o src/update_grape.o src/update_loc_p.o src/writebh.o \
 src/writerestartfile.o src/predict_potential.o src/wtimef.o 

GLOBJECTS=src/fthread.o src/viewer.o src/start_viewer.o src/libpt.a

ifndef CUDA_LIBDIRS
	CUDA_LIBDIRS = -L/data1/pelupes/libraries/CUDA/cuda/lib 
endif

SAPPORO_LIBDIRS ?= -L/data1/pelupes/libraries/sapporo_v1.5
GRAPE6_LIBDIRS ?= -L$(GRAPE6)/lib

LLIBSDIR = $(SAPPORO_LIBDIRS) $(CUDA_LIBDIRS) $(GRAPE6_LIBDIRS)
GPULIBS  = $(CUDA_LIBS) -lsapporo -L/opt/local/lib -lboost_thread-mt -lpthread -L/usr/lib -lstdc++.6
GRAPE6LIBS  ?= -lg6lx
LIBS     =  -lm -lpthread

LIBRARIES= $(LLIBSDIR) $(LIBS)

GLINCLUDE = -I/usr/include/GL \
-I/home/inti/libraries/f90gl-1.2.11/include/GL/

GLLIB = -L/home/inti/libraries/f90gl-1.2.11/lib \
-lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut  -lguide 

X11LIB ?= -L/usr/X11R6/lib64 -lXaw -lXt -lXmu -lXi -lXext -lX11

THREADLIB ?= -L. -lpthread src/libpt.a

PGLIBS ?= -L. -lpg5

G6LIBS ?= -L../../../../lib/g6 -lg6

CFLAGS += -O2 $(OPTS)
FFLAGS += -O2 $(OPTS)
FFLAGS += -I./src/

all: worker_code 

worker_code: worker_code.f90 phiGRAPE_code.o
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $(FFLAGS) $^ -o $@ $(OBJECTS) $(G6LIBS) $(LIBRARIES)

glworker_code: glworker_code.f90 phiGRAPE_code.o
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	make -C src gl "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $(FFLAGS) $^ -o $@ $(OBJECTS) $(G6OBJS) $(LIBRARIES) \
	$(GLOBJECTS) $(GLLIB)  $(X11LIB) $(THREADLIB) $(LIBS)

worker_code_gpu: worker_code.f90 phiGRAPE_code.o
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $(FFLAGS) $^ -o $@ $(OBJECTS) $(LIBRARIES) $(GPULIBS)
	
worker_code_grape: worker_code.f90 phiGRAPE_code.o
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $(FFLAGS) $^ -o $@ $(OBJECTS) $(LIBRARIES) $(GRAPE6LIBS)
	
worker_code_phantom_grape: worker_code.f90 phiGRAPE_code.o
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $(FFLAGS) $^ -o $@ $(OBJECTS) $(LIBRARIES) $(PGLIBS)

worker_code.f90: interface.py __init__.py
	$(CODE_GENERATOR) --type=f90 interface.py PhiGRAPEInterface -o $@

glworker_code.f90: interface.py __init__.py
	$(CODE_GENERATOR) --type=f90 interface.py PhiGRAPEInterfaceGL -o $@

clean:
	make -C src clean
	rm -f __init__.py *.pyc
	rm -f worker_code_gpu worker_code_grape
	rm -f worker_code worker_code.f90
	rm -f glworker_code glworker_code.f90

__init__.py:
	touch $@

.f.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

.F.o:
	$(FC) $(FFLAGS) -c $<
