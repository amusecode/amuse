OBJECTS = compact_chem_optimized.o
FC ?= gfortran
FFLAGS ?= -O2

LIBFILES= $(OBJECTS) hash.o amuse_helpers.o

# krome makefile has changed to be ifort, with targets for gfortran
#~ KROME_OPTIONS=-unsafe -iRHS -n networks/react_COthin -compiler $(FC) -useN -clean
KROME_OPTIONS=-unsafe -iRHS -n networks/react_COthin -useN -clean -noVerbatimFile
#  -compact

ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
AR=xiar cru
else
AR=ar crs
endif

ifeq ($(findstring gfortran, $(notdir $(FORTRAN))), gfortran)
FFLAGS += -ffree-line-length-none
endif

amuse_interface:  libchem.a

test:

# note the hard change of fortran compiler and options
krome/build/krome:
	-(cd krome; ./krome $(KROME_OPTIONS) )
	cd ..
	sed -i 's/^fc *=.*/fc = $(FC)/' krome/build/Makefile
	sed -i 's_^switch *= .*_switch = $(FFLAGS)_' krome/build/Makefile
	make -C krome/build/

libchem.a: krome/build/krome hash.o amuse_helpers.o
	$(AR) $@ $(filter-out krome/build/test.o, $(wildcard krome/build/*.o)) hash.o amuse_helpers.o

%.o: %.f90
	$(FC) $(FFLAGS) -Ikrome/build/ -c -o $@ $<

clean:
	rm -f *.o *.mod
	rm -f libchem.a
	rm -f krome/build/*

distclean: clean
	rm -f libchem.a
