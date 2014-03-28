# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk

##########################################################
CXX ?= mpicxx
CC ?= mpicc

MPICXX ?= mpicxx

CXXFLAGS ?= -Wall -g -O2
CFLAGS ?= -Wall -g -O2

LIBS = 
##########################################################
OBJS = main.o Sakura.o Timestep.o Interaction.o Particles.o Particle.o Twobody.o Diagnostics.o Communicator.o
EXEC = main.exe
CODELIB = libsakura.a
all: $(OBJS) $(EXEC) $(CODELIB)

##########################################################
main.exe: $(OBJS)
	$(MPICXX) $(CXXFLAGS) -o $(EXEC) $(OBJS) $(LIBS)

libsakura.a: $(OBJS)
	rm -f $@
	ar crs $@ $(OBJS)
	ranlib $@
##########################################################
main.o: main.cpp Sakura.o Timestep.o Interaction.o Particles.o Particle.o Twobody.o Communicator.o
	$(MPICXX) $(CXXFLAGS) -c main.cpp
Sakura.o: Sakura.h Sakura.cpp Timestep.o Interaction.o Particles.o Particle.o Twobody.o Communicator.o
	$(MPICXX) $(CXXFLAGS) -c Sakura.cpp
Twobody.o: Twobody.h Twobody.cpp
	$(CXX) $(CXXFLAGS) -c Twobody.cpp
Timestep.o: Timestep.h Timestep.cpp Particles.o Particle.o
	$(CXX) $(CXXFLAGS) -c Timestep.cpp
Interaction.o: Interaction.h Interaction.cpp Particles.o Particle.o Communicator.o
	$(MPICXX) $(CXXFLAGS) -c Interaction.cpp
Particles.o: Particles.h Particles.cpp Particle.o
	$(CXX) $(CXXFLAGS) -c Particles.cpp
Particle.o: Particle.h Particle.cpp
	$(CXX) $(CXXFLAGS) -c Particle.cpp
Diagnostics.o: Diagnostics.h Diagnostics.cpp
	$(CXX) $(CXXFLAGS) -c Diagnostics.cpp
Communicator.o: Communicator.h Communicator.cpp
	$(MPICXX) $(CXXFLAGS) -c Communicator.cpp
##########################################################
clean:
	rm -f *~
	rm -f *.o
	rm -f *.exe $(CODELIB)



