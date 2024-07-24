OBJECTS = compact_chem_optimized.o
FC ?= gfortran
FFLAGS ?= -O2

LIBFILES= $(OBJECTS) hash.o amuse_helpers.o

# krome makefile has changed to be ifort, with targets for gfortran
#~ KROME_OPTIONS=-unsafe -iRHS -n networks/react_COthin -compiler $(FC) -useN -clean
KROME_OPTIONS=-unsafe -iRHS -n networks/react_COthin -useN -clean -noVerbatimFile
#  -compact

FFLAGS += -ffree-line-length-none -fallow-argument-mismatch

amuse_interface:  libchem.a

test:

# note the hard change of fortran compiler and options
krome/build/krome:
	cd krome && ./krome $(KROME_OPTIONS)
	$(MAKE) -C krome/build/ fc="$(FC)" switch="$(FFLAGS)"

libchem.a: krome/build/krome hash.o amuse_helpers.o
	$(AR) crs $@ $(filter-out krome/build/test.o, $(wildcard krome/build/*.o)) hash.o amuse_helpers.o

%.o: %.f90 krome/build/krome
	$(FC) $(FFLAGS) -Ikrome/build/ -c -o $@ $<

clean:
	rm -f *.o *.mod
	rm -f libchem.a
	rm -f krome/build/*

distclean: clean
	rm -f libchem.a
