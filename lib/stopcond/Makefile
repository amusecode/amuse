AMUSE_DIR?=../..
include ${AMUSE_DIR}/config.mk

MPICC ?= mpicc

LIBNAME = libstopcond.a
LIBNAMEMPI = libstopcondmpi.a

OBJS = stopcond.o
OBJSF = stopcondf_isoc.o
OBJSMPI = stopcond.mo

all:$(LIBNAME) $(LIBNAMEMPI)

CFLAGS += -O2

$(LIBNAME): $(OBJS) $(OBJSF) 
	ar -r  $(LIBNAME) $^
	ranlib $(LIBNAME)
	
$(LIBNAMEMPI): $(OBJSMPI) $(OBJSF)
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
