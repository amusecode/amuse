# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPIFC ?= mpif90
FORTRAN ?= gfortran
FC ?= $(FORTRAN)


CODE_GENERATOR = $(AMUSE_DIR)/build.py

SRCDIR=src/trunk
BUILDDIR=build

CMAKE_DEPS=$(SRCDIR)/CMakeLists.txt $(SRCDIR)/cmake/CompilerFlags.cmake $(SRCDIR)/cmake/FileList.cmake

SRC = $(SRCDIR)/code/libint.f90

EVTWINLIB=$(SRCDIR)/lib/libtwin.a $(SRCDIR)/lib/libtwincore.a

FLIB= -L/software/local/intel64/compiler91/lib -lifcore -lifport
LD = ld
AR = ar ruv
RANLIB = ranlib
CMAKE_GENERATOR ?= Unix Makefiles

all:    $(BUILDDIR)/Makefile evtwin_worker

$(BUILDDIR)/Makefile: $(CMAKE_DEPS)
	-mkdir -p $(BUILDDIR)
	cd $(BUILDDIR) && FC="$(FC) $(FCFLAGS)" cmake ../src/trunk -DCMAKE_Fortran_FLAGS:STRING="$(FCFLAGS)"  -G "$(CMAKE_GENERATOR)" -DWANT_MKMERGERMOD=0 -DWANT_LIBRARY=1 

clean:
	rm -f *.pyc worker_code.f90 evtwin_worker
	make -C $(BUILDDIR) clean
	rm -rf $(BUILDDIR)
	rm -f $(EVTWINLIB)
	rm -rf $(SRCDIR)/lib

evtwin_worker: worker_code.f90 $(SRC)
	make -s --no-print-directory -C $(BUILDDIR)
	$(MPIFC) -I$(BUILDDIR) $(FS_FLAGS) $(LDFLAGS) $< -o $@ $(EVTWINLIB) $(FS_LIBS)  $(LIBS)

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py EVtwinInterface -o $@
