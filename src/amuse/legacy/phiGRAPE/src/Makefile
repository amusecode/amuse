
OPTS ?= -DNOMPI

CODELIB = libphigrape.a

AR = ar ruv
RANLIB = ranlib
RM = rm

AMUSE_DIR?=../../../../../

SOURCEF = \
allgather_ap.F      cmcorr.F        corrector.F     cputime.F		\
energy.F       	    get_min_t.F     gravity.F       initgrape.F		\
initmpi.F           initpars.F      output.F				\
predictor.F         readbodydata.F  readpars.F      readrestartfile.F	\
selectactive.F      selectsync.F    sendbodies2grape.F		    	\
sendbodydata.F      							\
sendinpdata.F       sumforces.F	    sync.F				\
timestep.F          update_grape.F  update_loc_p.F  writebh.F		\
writerestartfile.F predict_potential.F

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

SOURCEC = wtimef.c

SAF = main.F

CUDA_LIBDIRS ?= -L/data1/pelupes/programs/cuda/lib 
CUDA_LIBS ?= -lcudart
SAPPORO_LIBDIRS ?= -L/data1/pelupes/libraries/sapporo_v1.5


LLIBSDIR = $(SAPPORO_LIBDIRS) $(CUDA_LIBDIRS)
LIBS     =  $(CUDA_LIBS) -lsapporo -lm -lboost_thread-mt -lpthread


G6LIBS ?= -L../../../../../lib/g6 -lg6

LIBS     = -lm

LIBRARIES= $(LLIBSDIR) $(LIBS) $(G6LIBS) $(SCLIBS)

GLINCLUDE ?= -I/usr/include/GL \
-I/data/pelupes/libraries/f90gl-1.2.11/include/GL/

GLLIB = -L/data/pelupes/libraries/f90gl-1.2.11/lib \
-lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut 

X11LIB ?= -L/usr/X11R6/lib64 -lXaw -lXt -lXmu -lXi -lXext -lX11


THREAD_OPTIONS ?= -L. -lpthread -lpt

#------------------------------------------------------------
#---   C-definitions
#------------------------------------------------------------

CFLAGS  += -O3 $(OPTS) 

#------------------------------------------------------------
#---   FORTRAN-definitions
#------------------------------------------------------------

FFLAGS  += -O2 $(OPTS) 
#gfortran :-fno-second-underscore -fPIC -Wall

OBJECTS = $(SOURCEF:.F=.o) $(SOURCE:.f=.o) $(SOURCEC:.c=.o)
CODEOBJS = $(OBJECTS)
SAO     = $(SAF:.F=.o)

all: debuginfo phiGRAPE $(CODELIB)

debuginfo:
	@echo "enetered phigrape source makefile"
	@echo "OBJECTS  = " $(OBJECTS)
	@echo "CODEOBJS = " $(CODEOBJS)
	@echo $(SCLIBS)

gl: start_viewer.o viewer.o libpt.a fthread.o


$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

stars3d.o: stars3d.f90
	$(FC) -fpp $(FFLAGS) $(GLINCLUDE) -c -o $@ $<

viewer.o: viewer.f90
	$(FC) -fpp $(FFLAGS) $(GLINCLUDE) -c -o $@ $<

libpt.a: pt.c ptf77.c pt.h
	cc -c $(PIC) pt.c
	cc -c $(PIC) ptf77.c
	ar crv libpt.a pt.o ptf77.o
	ranlib libpt.a

phiGRAPE: $(OBJECTS) $(SAO) 
	$(FC) $(FFLAGS) -o phiGRAPE $^ $(LIBRARIES) 


clean:
	rm -f *.o .m* .n* *.lst *.L *~ phiGRAPE libpt.a *.mod $(CODELIB)

.f.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

.F.o:
	$(FC) $(FFLAGS) -c $<
