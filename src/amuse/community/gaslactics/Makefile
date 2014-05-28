# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

CC ?= gcc
CFLAGS = -O3

FC ?= gfortran
ifeq ($(findstring gfortran, $(shell $(FC) --help)), gfortran)
FFLAGS = -ffixed-line-length-0 -O0 -fno-backslash
FFLAGS2 = -O0
LIBS = -ldl -lpthread -lm -lstdc++
else
FFLAGS = -O0 -132 -nofor-main
FFLAGS2 = -O0 -nofor-main
LIBS = -L/opt/intel/Compiler/11.1/064/bin/intel64/ifort \
 -lifport -lifcore -limf -lsvml -lipgo -lintlc -ldl \
 -lpthread -lcxaguard -lm -lstdc++
endif

BINDIR = src/bin/
EXES = $(BINDIR)dbh $(BINDIR)genhalo $(BINDIR)genbulge \
  $(BINDIR)getfreqs $(BINDIR)gendisk $(BINDIR)diskdf \
  $(BINDIR)vcirc $(BINDIR)plothalo $(BINDIR)plotforce \
  $(BINDIR)gengas $(BINDIR)toascii


all: $(BINDIR) $(EXES) 

$(BINDIR):
	cd src; mkdir bin

$(EXES): Makefile
	make -C src/src all CC="$(FC)" CFLAGS="$(CFLAGS)" FC="$(FC)" FFLAGS="$(FFLAGS)" FFLAGS2="$(FFLAGS2)" LIBS="$(LIBS)"
	make -C src/src install

clean:
	rm -f *.pyc
	rm -f $(BINDIR)*
	make -C src/src purge
