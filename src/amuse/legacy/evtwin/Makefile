
FORTRAN ?= gfortran
AMUSE_DIR ?= ../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

LIBDIR=src/lib
OBJDIR=src/obj

LIBTWINA=$(LIBDIR)/libtwin.a
OBJ=$(OBJDIR)/twinmuse.o 
OBJ_EVT=obj/twinmuse.o 
MODDIR=src/modules
EVTWINLIB=src/lib/libtwin.a
 
FLIB= -L/software/local/intel64/compiler91/lib -lifcore -lifport

all:    src/makefile.dep src/makefile.sys worker_code

src/makefile.dep: src/Makefile
	@echo building makefile.deb
	@cd src && perl depend.pl code/*.f > makefile.dep
	@echo done
	

src/makefile.sys: src/Makefile
	@echo running configure
	@cd src && perl configure.pl > configure.out
	@echo done	
	
	
clean:
	make -s -C src clean
	/bin/rm -f *.pyc make.log
	/bin/rm -f mpi_code.f90 worker_code twin_library_v2.mod mesh.mod
	/bin/rm -f src/makefile.dep src/makefile.sys

distclean:
	make -C src distclean
	/bin/rm -f *.pyc

worker_code: mpi_code.f90 src/makefile.sys src/makefile.dep
	@-mkdir -p src/modules
	@-mkdir -p src/obj
	@-mkdir -p src/lib
	make -s --no-print-directory -C src lib/libtwin.a FORT=$(FORTRAN) >& make.log
	@cp $(MODDIR)/twin_library_v2.mod .
	@cp $(MODDIR)/mesh.mod .
	mpif90 $< -o $@ $(EVTWINLIB)

mpi_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py EVtwinInterface -o $@
	
.f.o: $<
	$(FORTRAN) $(FFLAGS) -c -o $@ $< 
