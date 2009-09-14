SOURCEF = \
allgather_ap.F      cmcorr.F        corrector.F     cputime.F		\
energy.F       	    get_min_t.F     gravity.F       initgrape.F		\
initmpi.F           initpars.F      output.F				\
predictor.F         readbodydata.F  readpars.F      readrestartfile.F	\
selectactive.F      selectsync.F    sendbodies2grape.F		    	\
sendbodydata.F      							\
sendinpdata.F       sumforces.F	    sync.F				\
timestep.F          update_grape.F  update_loc_p.F  writebh.F		\
writerestartfile.F

SOURCEC = wtimef.c g6lib.c 

SAF = main.F

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
LIBS     =  $(CUDA_LIBS) -lsapporo -lm -lboost_thread-mt -lpthread
LIBS     = -lm

LIBRARIES= $(LLIBSDIR) $(LIBS)

GLINCLUDE = -I/usr/include/GL \
-I/data1/pelupes/libraries/f90gl-1.2.11/include/GL/

GLLIB = -L/data1/pelupes/libraries/f90gl-1.2.11/lib \
-lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut 
GLLIB = 

X11LIB = -L/usr/X11R6/lib64 -lXaw -lXt -lXmu -lXi -lXext -lX11

THREAD_OPTIONS = -L. -lpthread -lpt

#------------------------------------------------------------
#---   C-definitions
#------------------------------------------------------------

# MPI
CC      = mpicc
CFLAGS  = -O2
# no MPI
#CC	= gcc
#CFLAGS  = -O2 -DNOMPI -fPIC
#PIC= -fPIC

#------------------------------------------------------------
#---   FORTRAN-definitions
#------------------------------------------------------------

# MPI
FORTRAN = mpif90
FFLAGS  = -O2 

# no MPI

#FC = ifort
#FFLAGS  = -O2 -DNOMPI -fPIC
#gfortran
#-fno-second-underscore -fPIC -Wall

OBJECTS = $(SOURCEF:.F=.o) $(SOURCE:.f=.o) $(SOURCEC:.c=.o)
SAO     = $(SAF:.F=.o)

all: phiGRAPE stars3d.o libpt.a fthread.o

stars3d.o: stars3d.f90
	$(FORTRAN) -fpp $(FFLAGS) $(GLINCLUDE) -c -o $@ $<

libpt.a: pt.c ptf77.c pt.h
	cc -c $(PIC) pt.c
	cc -c $(PIC) ptf77.c
	ar crv libpt.a pt.o ptf77.o
	ranlib libpt.a

phiGRAPE: $(OBJECTS) $(SAO) 
	$(FORTRAN) $(FFLAGS) -o phiGRAPE $^ $(LIBRARIES) 


clean:
	rm -f *.o .m* .n* *.lst *.L *~ phiGRAPE

.f.o:
	$(FORTRAN) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

.F.o:
	$(FORTRAN) $(FFLAGS) -c $<
