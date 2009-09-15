
FC = mpif90
CC = mpicc

OPTS+=-DNOMPI

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
 src/writerestartfile.o src/wtimef.o src/g6lib.o

GLOBJECTS=src/fthread.o src/stars3d.o src/libpt.a

ifndef CUDA_LIBDIRS
	CUDA_LIBDIRS = -L/data1/pelupes/programs/cuda/lib 
endif

ifndef CUDA_LIBS
	CUDA_LIBS = -lcudart
endif


ifndef SAPPORO_LIBDIRS
	SAPPORO_LIBDIRS = -L/data1/pelupes/libraries/sapporo_v1.5
endif


LLIBSDIR = $(SAPPORO_LIBDIRS) $(CUDA_LIBDIRS)
LIBS     = $(CUDA_LIBS) -lsapporo -lboost_thread-mt -lpthread 
LIBS     = -lm  
LIBS     =  -lm -lpthread

LIBRARIES= $(LLIBSDIR) $(LIBS)

GLINCLUDE = -I/usr/include/GL \
-I/data1/pelupes/libraries/f90gl-1.2.11/include/GL/

GLLIB = -L/data1/pelupes/libraries/f90gl-1.2.11/lib \
-lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut  -lguide 

X11LIB = -L/usr/X11R6/lib64 -lXaw -lXt -lXmu -lXi -lXext -lX11

THREADLIB = -L. -lpthread src/libpt.a


CFLAGS  = -O2 $(OPTS) -fPIC

FFLAGS  = -O2 $(OPTS) -fPIC
# gfortran: -fno-second-underscore -fPIC -Wall

FFLAGS+= -I./src/


all: muse_worker


muse_worker: muse_worker.f90 muse_interface.o 
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $^ -o $@ $(OBJECTS) $(LIBRARIES)

muse_worker_gl: muse_worker_gl.f90 muse_interface.o muse_interface_gl.o 
	make -C src all "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	make -C src gl "FC=$(FC)" "OPTS=$(OPTS)" "CC=$(CC)"
	mpif90 $^ -o $@ $(OBJECTS) $(GLOBJECTS) $(LIBRARIES) $(GLLIB) $(X11LIB) $(THREADLIB)
        
muse_worker.f90: muse_dynamics_mpi.py
	../../../../bin/create_f90_worker.py muse_dynamics_mpi.py PhiGRAPE > $@

muse_worker_gl.f90: muse_dynamics_gl.py
	../../../../bin/create_f90_worker.py muse_dynamics_gl.py PhiGRAPE > $@


interface_gl:	$(OBJECTS) $(MUSE_IF) muse_interface_gl.F 


clean:
	make -C src clean
	rm -f *.o .m* .n* *.lst *.L *~ muse_worker
	rm muse_worker muse_worker.f90 muse_worker_gl.f90

.f.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

.F.o:
	$(FC) $(FFLAGS) -c $<
