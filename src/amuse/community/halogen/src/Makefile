# Makefile for halogen4muse

# Library for AMUSE

CODELIB = libhalogen.a

# Executable

BASE    = halogen4muse
EXT     = 
EXE     = $(BASE)$(EXT)

# Compiler stuff

CC	?= gcc
CFLAGS	= -O3 -Wall
LIBS	= -lm

# Object definition

LIBOBJS = functions.o amuse_routines.o
EXEOBJS = $(BASE).o functions.o routines.o

# Rules

all: $(EXE) $(CODELIB)

$(EXE):	$(EXEOBJS) Makefile
	$(CC) $(CFLAGS) $(EXEOBJS) -o $(EXE) $(LIBS)

$(CODELIB): $(LIBOBJS) Makefile
	ar crs $@ $(LIBOBJS)

clean:
	-rm -f *.o *~ $(EXE) $(CODELIB)

tar:
	cd ..; tar cvf - $(BASE)/*.c $(BASE)/*.h $(BASE)/Makefile > $(BASE).tar

# Dependencies

halogen4muse.o: definitions.h functions.h routines.h
functions.o: definitions.h functions.h routines.h
routines.o: definitions.h functions.h routines.h

amuse_routines.o: routines.c definitions.h functions.h routines.h
	$(CC) $(CFLAGS) -DAMUSE_FLAG -c -o $@ $<
