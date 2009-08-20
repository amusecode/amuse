
CC = gcc
CPP = g++

CFLAGS = -g 
CPPFLAGS = -Wall -g -DTREE -DNOGRAPHICS 
#CPPFLAGS = -Wall -g -DNOGRAPHICS

OBJ = BHTC.o nbody.o BHtree.o gravity.o pgetopt.o second.o


all: BHTC

clean:
	rm -f *.o BHTC *~


BHTC: $(OBJ)
	$(CPP) $(CPPFLAGS) $(OBJ) -lm -o BHTC

.c.o:
	$(CC) $(CFLAGS) -c $?

.C.o:
	$(CPP) $(CPPFLAGS) -c $?

test: compile

compile: BHTC
