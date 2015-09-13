AMUSE_DIR?=../..
include ${AMUSE_DIR}/config.mk

MPICC ?= mpicc

LIBNAME = libstopcond.a
LIBNAMEMPI = libstopcondmpi.a

OBJS = stopcond.o
OBJSMPI = stopcond.mo
ifeq ($(FC_ISO_C_AVAILABLE), yes)
OBJS += stopcondf_isoc.o
OBJSMPI += stopcondf_isoc.o
endif

all:$(LIBNAME) $(LIBNAMEMPI)

CFLAGS += -O2

$(LIBNAME): $(OBJS)
	ar -r  $(LIBNAME) $^
	ranlib $(LIBNAME)
	
$(LIBNAMEMPI): $(OBJSMPI)
	ar -r  $(LIBNAMEMPI) $^
	ranlib $(LIBNAMEMPI)

%.mo: %.c
	$(MPICC) $(CFLAGS) -DMPILIB  -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
	
%.o: %.F90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -Rf *.o *.lo *.a *.mo *.mod
