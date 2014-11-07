OBJECTS = compact_chem_optimized.o
FC ?= gfortran
FFLAGS ?= -O3

LIBFILES= $(OBJECTS) hash.o amuse_helpers.o

KROME_OPTIONS=-unsafe -iRHS -n networks/react_COthin -compiler $(FC) -useN -clean
#  -compact

ifeq ($(FC), ifort)
AR=xiar cru
else
AR=ar crs
endif

amuse_interface:  libchem.a

test:

krome/build/krome:
	-(cd krome; ./krome $(KROME_OPTIONS) )
	cd ..
	make -C krome/build/

libchem.a: krome/build/krome hash.o amuse_helpers.o
	$(AR) $@ krome/build/[!test]*.o hash.o amuse_helpers.o

%.o: %.f90
	$(FC) $(FFLAGS) -Ikrome/build/ -c -o $@ $<

clean:
	rm -f *.o *.mod
	rm -f libchem.a
	rm -f krome/build/*

distclean: clean
	rm -f libchem.a
