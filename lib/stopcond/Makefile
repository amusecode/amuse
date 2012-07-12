AMUSE_DIR?=../..
include ${AMUSE_DIR}/config.mk

MPICC ?= mpicc

LIBNAME = libstopcond.a
LIBFNAME = libstopcondf.a
LIBNAMEMPI = libstopcondmpi.a
LIBFNAMEMPI = libstopcondfmpi.a

OBJS = stopcond.o
OBJSF = stopcondf.o
OBJSMPI = stopcond.mo

all:$(LIBNAME) $(LIBNAMEMPI) $(LIBFNAME)

CFLAGS += -O2

$(LIBNAME): $(OBJS)
	ar -r  $(LIBNAME) $(OBJS)
	ranlib $(LIBNAME)
	
$(LIBFNAME): $(OBJSF)
	ar -r  $(LIBFNAME) $(OBJSF)
	ranlib $(LIBFNAME)

$(LIBNAMEMPI): $(OBJSMPI)
	ar -r  $(LIBNAMEMPI) $(OBJSMPI)
	ranlib $(LIBNAMEMPI)


%.mo: %.c
	$(MPICC) $(CFLAGS) -DMPILIB  -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
	
%.o: %.F90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -Rf *.o *.lo *.a *.mo
