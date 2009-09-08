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
 src/writerestartfile.o src/wtimef.o

LLIBSDIR = -L/data1/pelupes/libraries/sapporo_v1.5 \
-L/data1/pelupes/programs/cuda/lib 
LIBS     = -lguide -lsapporo -lm -lcudart -lboost_thread-mt

LIBRARIES= $(LLIBSDIR) $(LIBS)

GLINCLUDE = -I/usr/include/GL \
-I/data1/pelupes/libraries/f90gl-1.2.11/include/GL/

GLLIB = -L/data1/pelupes/libraries/f90gl-1.2.11/lib \
-lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut 

X11LIB = -L/usr/X11R6/lib64 -lXaw -lXt -lXmu -lXi -lXext -lX11

THREAD_OPTIONS = -L. -lpthread -lpt

CC	= gcc
CFLAGS  = -O2 -DNOMPI -fPIC
PIC= -fPIC

# MPI
#FORTRAN = mpif77
#FFLAGS  = -O2 -Wno-globals
FORTRAN = ifort
FC = ifort
FFLAGS  = -O2 -DNOMPI -fPIC
# gfortran -fno-second-underscore -fPIC -Wall

FFLAGS+= -I./src/


all: muse_worker


muse_worker: muse_worker.f90 muse_interface.o 
	make -C src all
	mpif90 $^ -o $@ $(OBJECTS) $(LIBRARIES)
        
muse_worker.f90: muse_dynamics_mpi.py
	../../../../bin/create_f90_worker.py muse_dynamics_mpi.py PhiGRAPE > $@


interface_gl:	$(OBJECTS) $(MUSE_IF) muse_interface_gl.F fthread.o stars3d.o libpt.a


clean:
	make -C src clean
	rm -f *.o .m* .n* *.lst *.L *~

.f.o:
	$(FORTRAN) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

.F.o:
	$(FORTRAN) $(FFLAGS) -c $<
