
MPIF90 ?= mpif90
FORTRAN ?= gfortran
FC ?= $(FORTRAN)
CMAKE_Fortran_COMPILER ?= $(FORTRAN)

AMUSE_DIR ?= ../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

LIBDIR=src/lib
SRCDIR=src/code
MODDIR=src/modules

CMAKE_DEPS=src/CMakeLists.txt src/cmake/CompilerFlags.cmake src/cmake/FileList.cmake

SRC=$(SRCDIR)/library_v2.f90
EVTWINLIB=$(LIBDIR)/libtwin.a
 

FLIB= -L/software/local/intel64/compiler91/lib -lifcore -lifport
LD = ld
AR = ar ruv
RANLIB = ranlib

all:    create_Makefile evtwin_worker

create_Makefile: $(CMAKE_DEPS)
	mkdir -p $(MODDIR)
	cd $(MODDIR) && cmake ..

clean:
	-make -C $(MODDIR) clean
	rm -f *.pyc make.log
	rm -f mpi_code.f90 worker_code evtwin_worker twin_library_v2.mod mesh.mod
	rm -f $(MODDIR)/CMakeCache.txt
	rm -rf $(MODDIR)/CMakeFiles
	rm -rf $(MODDIR)/Makefile

distclean:
	make -C src distclean
	/bin/rm -f *.pyc

evtwin_worker: mpi_code.f90 $(SRC)
	@-mkdir -p $(LIBDIR)
#	make -s --no-print-directory -C src lib/libtwin.a FORT=$(FORTRAN)
	make -s --no-print-directory -C $(MODDIR)
#	cp $(MODDIR)/twin_library_v2.mod .
#	cp $(MODDIR)/mesh.mod .
	$(MPIF90) -I$(MODDIR) $< -o $@ $(EVTWINLIB)

mpi_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py EVtwinInterface -o $@
	
