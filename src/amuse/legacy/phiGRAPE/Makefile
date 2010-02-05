
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

G6OBJS = src/g6lib.o

GLOBJECTS=src/fthread.o src/viewer.o src/libpt.a

CUDA_LIBDIRS ?= -L/data1/pelupes/programs/cuda/lib 
CUDA_LIBS ?= -lcudart

SAPPORO_LIBDIRS ?= -L/data1/pelupes/libraries/sapporo_v1.5
GRAPE6_LIBDIRS ?= -L$(GRAPE6)/lib

LLIBSDIR = $(SAPPORO_LIBDIRS) $(CUDA_LIBDIRS) $(GRAPE6_LIBDIRS)
GPULIBS  = $(CUDA_LIBS) -lsapporo -L/opt/local/lib -lboost_thread-mt -lpthread -L/usr/lib -lstdc++.6
GRAPE6LIBS  ?= -lg6lx
LIBS     =  -lm -lpthread

LIBRARIES= $(LLIBSDIR) $(LIBS)

GLINCLUDE = -I/usr/include/GL \
-I/data1/pelupes/libraries/f90gl-1.2.11/include/GL/

GLLIB = -L/data1/pelupes/libraries/f90gl-1.2.11/lib \
-lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut  -lguide 

X11LIB ?= -L/usr/X11R6/lib64 -lXaw -lXt -lXmu -lXi -lXext -lX11

THREADLIB ?= -L. -lpthread src/libpt.a


CFLAGS += -O2 $(OPTS)
FFLAGS += -O2 $(OPTS)
FFLAGS += -I./src/


all: muse_worker worker_code 

muse_worker: muse_worker.f90 muse_interface.o
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $(FFLAGS) $^ -o $@ $(OBJECTS) $(G6OBJS) $(LIBRARIES)

worker_code: worker_code.f90 phiGRAPE_code.o
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $(FFLAGS) $^ -o $@ $(OBJECTS) $(G6OBJS) $(LIBRARIES)
	
worker_code_gpu: worker_code.f90 phiGRAPE_code.o
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $(FFLAGS) $^ -o $@ $(OBJECTS) $(LIBRARIES) $(GPULIBS)
	
worker_code_grape: worker_code.f90 phiGRAPE_code.o
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $(FFLAGS) $^ -o $@ $(OBJECTS) $(LIBRARIES) $(GRAPE6LIBS)

muse_worker_gpu: muse_worker.f90 muse_interface.o 
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $(FFLAGS) $^ -o $@ $(OBJECTS) $(LIBRARIES) $(GPULIBS)
	
muse_worker_grape: muse_worker.f90 muse_interface.o 
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $(FFLAGS) $^ -o $@ $(OBJECTS) $(LIBRARIES) $(GRAPE6LIBS)


muse_worker_gl: muse_worker_gl.f90 muse_interface.o muse_interface_gl.o 
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	make -C src gl "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $^ -o $@ $(OBJECTS) $(G6OBJS) $(GLOBJECTS) $(LIBRARIES) \
	  $(GLLIB)  $(X11LIB) $(THREADLIB)
        
muse_worker.f90: muse_dynamics_mpi.py
	$(CODE_GENERATOR) --type=f90 muse_dynamics_mpi.py PhiGRAPE -o $@

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py PhiGRAPE -o $@

muse_worker_gl.f90: muse_dynamics_gl.py
	$(CODE_GENERATOR) --type=f90 muse_dynamics_gl.py PhiGRAPE -o $@


interface_gl:	$(OBJECTS) $(MUSE_IF) muse_interface_gl.F 


clean:
	make -C src clean
	rm -f *.o .m* .n* *.lst *.L *~ muse_worker muse_worker.f90 
	rm -f muse_worker_gl muse_worker_gl.f90
	rm -f muse_worker_gpu muse_worker_grape 
	rm -f worker_code_gpu worker_code_grape
	rm -f worker_code worker_code.f90


.f.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

.F.o:
	$(FC) $(FFLAGS) -c $<
