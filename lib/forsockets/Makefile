LIBNAME = libforsockets.a
LIBFNAME = libforsocketsf.a
OBJS = forsockets.o
OBJSF = forsocketsf.o

MPIF90 ?= openmpif90
MPICC ?= openmpicc

all:$(LIBNAME) $(LIBFNAME)

CFLAGS += -O2

$(LIBNAME): $(OBJS)
	ar -r  $(LIBNAME) $(OBJS)
	ranlib $(LIBNAME)

$(LIBFNAME): $(OBJSF)
	ar -r  $(LIBFNAME) $(OBJSF)
	ranlib $(LIBFNAME)

%.o: %.c
	$(CC) $(CFLAGS) -Wall -c -o $@ $<
	
%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -Rf *.o *.lo *.a *.mo *.mod
