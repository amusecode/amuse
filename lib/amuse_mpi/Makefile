AMUSE_DIR?=../..
include ${AMUSE_DIR}/config.mk

MPICC ?= mpicc

LIBNAME = libamuse_mpi.a

OBJS = amuse_mpi.o
#ifeq ($(FC_ISO_C_AVAILABLE), yes)
#OBJS += stopcondf_isoc.o
#OBJSMPI += stopcondf_isoc.o
#endif

all:$(LIBNAME)

CFLAGS += -fPIC

$(LIBNAME): $(OBJS)
	ar -r  $(LIBNAME) $^
	ranlib $(LIBNAME)
	
%.o: %.c
	$(MPICC) $(CFLAGS) -c -o $@ $<
	
%.o: %.F90
	$(MPIFC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -Rf *.o *.a *.mod
