# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(strip $(AMUSE_DIR)),)
    AMUSE_DIR = ../../../..
endif
-include ${AMUSE_DIR}/config.mk

ifeq ($(strip $(CC)),)
    CC = gcc
endif
ifeq ($(strip $(CFLAGS)),)
    CFLAGS = -O3
endif

ifeq ($(strip $(FC)),)
    FC = gfortran
endif
ifeq ($(strip $(FORTRAN)),)
    FORTRAN = $(FC)
endif
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

GASBINDIR = gas_src/bin/
GASEXES = $(GASBINDIR)dbh $(GASBINDIR)genhalo $(GASBINDIR)genbulge  \
  $(GASBINDIR)getfreqs $(GASBINDIR)gendisk $(GASBINDIR)diskdf $(GASBINDIR)gengas $(GASBINDIR)toascii

all: $(BINDIR) $(EXES) $(GASBINDIR) $(GASEXES) 

$(BINDIR):
	cd src; mkdir bin

$(GASBINDIR):
	cd gas_src; mkdir bin

$(EXES): Makefile
	echo $(findstring xlf, $(notdir $(FORTRAN)))
	make -C src/src all CC="$(CC)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" F77="$(FC)" FFLAGS="$(FFLAGS)" F77FLAGS="$(F77FLAGS)"
	make -C src/src install
	make -C src/potsrc all CC="$(CC)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" F77="$(FC)" FFLAGS="$(FFLAGS)" F77FLAGS="$(F77FLAGS)"
	make -C src/potsrc install

$(GASEXES): Makefile 
	echo $(findstring xlf, $(notdir $(FORTRAN)))
	make -C gas_src/src all CC="$(CC)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" F77="$(FC)" FFLAGS="$(FFLAGS)" FFLAGS2="$(FCFLAGS)" F77FLAGS="$(F77FLAGS)" PGSTUB=1
	make -C gas_src/src install


clean:
	rm -f *.pyc
	make -C src/src clean
	make -C gas_src/src purge
	make -C src/potsrc clean
	rm -f src/bin/*
	rm -f gas_src/bin/*
