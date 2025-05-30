# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk

###################################################################

CXX ?= g++
CC ?= gcc

CXXFLAGS ?= -Wall -g -O2
CXXFLAGS += -std=c++11
CFLAGS ?= -Wall -g -O2

GMP_LIBS ?= -lgmp
MPFR_LIBS ?= -lmpfr
LIBS = $(MPFR_LIBS) $(GMP_LIBS)

###################################################################

OBJS = integrator_MS.o Bs_integrator.o Cluster.o Clock.o Dynamics.o Force.o Star.o  dlmalloc.o
EXEC = integrator_MS.exe
CODELIB = libadaptb.a
all: $(OBJS) $(EXEC) $(CODELIB)

###################################################################

integrator_MS.exe: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJS) $(LIBS)

libadaptb.a: $(OBJS)
	rm -f $@
	ar crs $@ $(OBJS)
	ranlib $@

###################################################################

integrator_MS.o: integrator_MS.cpp Bs_integrator.o
	$(CXX) $(CXXFLAGS) $(MPFR_FLAGS) -c integrator_MS.cpp

Bs_integrator.o: Bs_integrator.cpp Bs_integrator.h Cluster.o
	$(CXX) $(CXXFLAGS) $(MPFR_FLAGS) -c Bs_integrator.cpp

Cluster.o: Cluster.cpp Cluster.h Star.o Force.o Dynamics.o
	$(CXX) $(CXXFLAGS) $(MPFR_FLAGS) -c Cluster.cpp

Clock.o: Clock.cpp Clock.h dlmalloc.o
	$(CXX) $(CXXFLAGS) $(MPFR_FLAGS) -c Clock.cpp

Dynamics.o: Dynamics.cpp Dynamics.h Star.o
	$(CXX) $(CXXFLAGS) $(MPFR_FLAGS) -c Dynamics.cpp

Force.o: Force.cpp Force.h Star.o
	$(CXX) $(CXXFLAGS) $(MPFR_FLAGS) -c Force.cpp

Star.o: Star.cpp Star.h dlmalloc.o
	$(CXX) $(CXXFLAGS) $(MPFR_FLAGS) -c Star.cpp

dlmalloc.o: dlmalloc.c dlmalloc.h
	$(CC) $(CFLAGS) -c dlmalloc.c

###################################################################

clean:
	rm -f *~
	rm -f *.o
	rm -f *.exe $(CODELIB)



