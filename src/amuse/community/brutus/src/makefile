# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk
###################################################################
## Change these paths

#HOME    = /data1/boekholt/
#MPRDIR  = $(HOME)Packages/mpreal_local/include/ # mpreal.h
#MPRLIB  = $(HOME)Packages/mpreal_local/lib/	# libmpfr.a, libgmp.a

###################################################################
COMP = $(CXX)
CFLAGS = -O1 $(CXXFLAGS)
CC  = $(COMP) $(CFLAGS) -c
CO  = $(COMP) $(CFLAGS) -o
###################################################################
OBJS = main.o Brutus.o Bulirsch_Stoer.o Cluster.o Star.o
EXE = main.exe
CODELIB = libbrutus.a

all: $(OBJS) $(EXEC) $(CODELIB)
###################################################################
main.exe: main.o Brutus.o Bulirsch_Stoer.o Cluster.o Star.o  
	$(CO) main.exe main.o Brutus.o Bulirsch_Stoer.o Cluster.o Star.o $(MPRLIB)libmpfr.a $(MPRLIB)libgmp.a

libbrutus.a: $(OBJS)
	rm -f $@
	ar crs $@ $(OBJS)
	ranlib $@
###################################################################
main.o: main.cpp Brutus.o Bulirsch_Stoer.o Cluster.o Star.o 
	$(CC) main.cpp
###################################################################
Brutus.o: Brutus.h Brutus.cpp Bulirsch_Stoer.o Cluster.o Star.o 
	$(CC) Brutus.cpp

Bulirsch_Stoer.o: Bulirsch_Stoer.h Bulirsch_Stoer.cpp Cluster.o Star.o 
	$(CC) Bulirsch_Stoer.cpp 
Cluster.o: Cluster.h Cluster.cpp Star.o 
	$(CC) Cluster.cpp 
Star.o: Star.h Star.cpp 
	$(CC) Star.cpp 
###################################################################
clean:
	rm -f *~
	rm -f *.o
	rm -f *.exe $(CODELIB)
	rm -f -r ccache
	rm -f .fuse*



