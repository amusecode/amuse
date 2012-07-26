OBJECTS = makefractal.o
FC ?= gfortran
FFLAGS ?= -O3

LIBFILES= $(OBJECTS) amuse_helpers.o

amuse_interface:  libfractal.a

libfractal.a: $(LIBFILES)
	ar crs $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.mod

purge: clean
