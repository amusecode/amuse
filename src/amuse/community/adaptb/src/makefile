# path to libmpfr.a and libgmp.a (platform dependant)
LDIR = /home/boekholt/Packages/
# path to dlmalloc.o, dlmalloc.c, mpreal.cpp, mpreal.h
LLDIR = /home/boekholt/Packages/mpfr_c++/

###################################################################

CC  = g++ -O4 -c
CO  = g++ -O4 -o

###################################################################

OBJS = integrator_MS.o Bs_integrator.o Cluster.o Clock.o Dynamics.o Force.o Star.o mpreal.o dlmalloc.o
EXEC=integrator_MS.exe
CODELIB=libadaptb.a
all: $(OBJS) $(EXEC) $(CODELIB)

###################################################################

integrator_MS.exe: integrator_MS.o Bs_integrator.o Cluster.o Clock.o Dynamics.o Force.o Star.o mpreal.o dlmalloc.o
	$(CO) integrator_MS.exe integrator_MS.o Bs_integrator.o Cluster.o Clock.o Dynamics.o Force.o Star.o mpreal.o dlmalloc.o $(LDIR)libmpfr.a $(LDIR)libgmp.a

libadaptb.a: $(OBJS)
	rm -f $@
	ar crs $@ $(OBJS)
	ranlib $@

###################################################################

integrator_MS.o: integrator_MS.cpp Bs_integrator.o
	$(CC) -I${LLDIR} integrator_MS.cpp

Bs_integrator.o: Bs_integrator.cpp Bs_integrator.h Cluster.o
	$(CC) -I${LLDIR} Bs_integrator.cpp

Cluster.o: Cluster.cpp Cluster.h Star.o Force.o Dynamics.o
	$(CC) -I${LLDIR} Cluster.cpp

Clock.o: Clock.cpp Clock.h mpreal.o dlmalloc.o
	$(CC) -I${LLDIR} Clock.cpp

Dynamics.o: Dynamics.cpp Dynamics.h Star.o
	$(CC) -I${LLDIR}  Dynamics.cpp

Force.o: Force.cpp Force.h Star.o
	$(CC) -I${LLDIR} Force.cpp

Star.o: Star.cpp Star.h mpreal.o dlmalloc.o
	$(CC) -I${LLDIR} Star.cpp

mpreal.o: $(LLDIR)mpreal.cpp $(LLDIR)mpreal.h
	$(CC) -I${LLDIR} $(LLDIR)mpreal.cpp

dlmalloc.o: $(LLDIR)dlmalloc.c $(LLDIR)dlmalloc.h
	$(CC) -I${LLDIR}  $(LLDIR)dlmalloc.c

###################################################################

clean:
	rm -f *~
	rm -f *.o
	rm -f *.exe $(CODELIB)



