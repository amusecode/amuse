# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

CC ?= gcc
CFLAGS ?= -O3

FC ?= gfortran
FORTRAN ?= $(FC)
ifeq ($(findstring gfortran, $(shell $(FORTRAN) --help)), gfortran)
FFLAGS = $(FCFLAGS) -ffixed-line-length-0 -O -fno-backslash
endif
ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
FFLAGS = $(FCFLAGS) -O3 -132 -nofor-main
endif

ifneq (,$(findstring xlf, $(notdir $(FORTRAN))))
FFLAGS = $(FCFLAGS) -qfixed -qextname
endif


BINDIR = src/bin/
EXES = $(BINDIR)dbh $(BINDIR)genhalo $(BINDIR)genbulge \
  $(BINDIR)getfreqs $(BINDIR)gendisk $(BINDIR)diskdf $(BINDIR)toascii


all: $(BINDIR) $(EXES) 

$(BINDIR):
	cd src; mkdir bin

$(EXES): Makefile
	echo $(findstring xlf, $(notdir $(FORTRAN)))
	make -C src/src all CC="$(CC)" CFLAGS="$(CFLAGS)" F77="$(FC)" FFLAGS="$(FFLAGS)" F77FLAGS="$(F77FLAGS)"
	make -C src/src install
	make -C src/potsrc all CC="$(CC)" CFLAGS="$(CFLAGS)" F77="$(FC)" FFLAGS="$(FFLAGS)" F77FLAGS="$(F77FLAGS)"
	make -C src/potsrc install

clean:
	rm -f *.pyc
	make -C src/src clean
	make -C src/potsrc clean

