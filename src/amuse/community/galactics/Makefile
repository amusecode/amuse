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


BINDIR = data/bin/
EXES = $(BINDIR)dbh $(BINDIR)genhalo $(BINDIR)genbulge \
  $(BINDIR)getfreqs $(BINDIR)gendisk $(BINDIR)diskdf $(BINDIR)toascii

GASBINDIR = data/gbin/
GASEXES = $(GASBINDIR)dbh $(GASBINDIR)genhalo $(GASBINDIR)genbulge  \
  $(GASBINDIR)getfreqs $(GASBINDIR)gendisk $(GASBINDIR)diskdf $(GASBINDIR)gengas $(GASBINDIR)toascii

all: $(BINDIR) $(EXES) $(GASBINDIR) $(GASEXES) 

$(BINDIR):
	mkdir -p $(BINDIR)

$(GASBINDIR):
	mkdir -p $(GASBINDIR)

$(EXES): Makefile
	mkdir -p src/bin
	mkdir -p gas_src/bin
	echo $(findstring xlf, $(notdir $(FORTRAN)))
	make -C src/src all CC="$(CC)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" F77="$(FC)" FFLAGS="$(FFLAGS)" F77FLAGS="$(F77FLAGS)"
	make -C src/src install
	make -C src/potsrc all CC="$(CC)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" F77="$(FC)" FFLAGS="$(FFLAGS)" F77FLAGS="$(F77FLAGS)"
	make -C src/potsrc install
	cp src/bin/* $(BINDIR)/


$(GASEXES): Makefile 
	echo $(findstring xlf, $(notdir $(FORTRAN)))
	make -C gas_src/src all CC="$(CC)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" F77="$(FC)" FFLAGS="$(FFLAGS)" FFLAGS2="$(FCFLAGS)" F77FLAGS="$(F77FLAGS)" PGSTUB=1
	make -C gas_src/src install
	cp gas_src/bin/* $(GASBINDIR)/


clean:
	make -C src/src clean
	make -C gas_src/src purge
	make -C src/potsrc clean
	rm -f src/bin/*
	rm -f gas_src/bin/*

distclean: clean
	rm -f *.pyc
	rm -rf data

