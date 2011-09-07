CC ?= gcc
CFLAGS = -O3

FC ?= gfortran
ifeq ($(findstring gfortran, $(shell $(FC) --help)), gfortran)
FFLAGS = -ffixed-line-length-0 -O -fno-backslash
else
FFLAGS = -O3 -132 -nofor-main
endif

BINDIR = src/bin/
EXES = $(BINDIR)dbh $(BINDIR)genhalo $(BINDIR)genbulge \
  $(BINDIR)getfreqs $(BINDIR)gendisk $(BINDIR)diskdf $(BINDIR)toascii


all: $(BINDIR) $(EXES) 

$(BINDIR):
	cd src; mkdir bin

$(EXES): Makefile
	make -C src/src all CC="$(CC)" CFLAGS="$(CFLAGS)" F77="$(FC)" FFLAGS="$(FFLAGS)"
	make -C src/src install
	make -C src/potsrc all CC="$(CC)" CFLAGS="$(CFLAGS)" F77="$(FC)" FFLAGS="$(FFLAGS)"
	make -C src/potsrc install

clean:
	rm -f *.pyc
	make -C src/src clean
	make -C src/potsrc clean
