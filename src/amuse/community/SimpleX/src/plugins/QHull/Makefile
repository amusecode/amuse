# Unix Makefile for qhull --- modified for SimpleX/AMUSE
#
#       see README.txt
#       the Unix distribution contains a configure Makefile
#
#       make           to produce qhull library libqhull.a
#       make clean     to remove object files and library
#
#       HFILES -- .h files for printing
#       OBJS -- specifies the object files of libqhull.a
#

MPICC ?= mpicc
CFLAGS ?= -g -O2 -Wall

# OBJS in execution frequency order.  CFILES after qhull.c are alphabetical
OBJS = user.o global.o stat.o io.o geom2.o poly2.o \
       merge.o qhull.o geom.o poly.o qset.o mem.o

HFILES= user.h qhull.h qhull_a.h geom.h io.h mem.h merge.h poly.h qset.h stat.h

all: libqhull.a

%.o : %.c $(HFILES)
	$(MPICC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f *.o libqhull.a

libqhull.a: $(OBJS)
	ar r libqhull.a $(OBJS)

