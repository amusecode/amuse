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
OBJSFMPI = stopcondf.mo

all:$(LIBNAME) $(LIBNAMEMPI) $(LIBFNAME) $(LIBFNAMEMPI)

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

$(LIBFNAMEMPI): $(OBJSFMPI)
	ar -r  $(LIBFNAMEMPI) $(OBJSFMPI)
	ranlib $(LIBFNAMEMPI)

%.mo: %.c
	$(MPICC) $(CFLAGS) -DMPILIB  -c -o $@ $<

%.mo: %.F90
	$(MPICC) $(FCFLAGS) -DMPILIB  -c -o $@ $<


%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
	
%.o: %.F90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -Rf *.o *.lo *.a *.mo
