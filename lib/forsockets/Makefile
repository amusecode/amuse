AMUSE_DIR?=../..
include ${AMUSE_DIR}/config.mk

LIBNAME = libforsockets.a
OBJS = forsockets.o
ifeq ($(FC_ISO_C_AVAILABLE), yes)
OBJS += forsocketsf.o
endif

MPIF90 ?= mpif90
MPICC ?= openmpicc

all:$(LIBNAME) 

CFLAGS += -O2

$(LIBNAME): $(OBJS)
	ar -r  $(LIBNAME) $^
	ranlib $(LIBNAME)

%.o: %.c
	$(CC) $(CFLAGS) -Wall -c -o $@ $<
	
%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -Rf *.o *.lo *.a *.mo *.mod
