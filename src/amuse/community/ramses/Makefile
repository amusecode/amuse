# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPIFC ?= mpif90
FC ?= gfortran
CODE_GENERATOR = $(AMUSE_DIR)/build.py

SRC=src/bin

NVECTOR = 64
NPRE = 8
NVAR = 3
SOLVER = hydro
DEFINES1D = -DNVECTOR=$(NVECTOR) -DNDIM=1 -DNPRE=$(NPRE) -DNVAR=$(NVAR) -DSOLVER$(SOLVER)
DEFINES2D = -DNVECTOR=$(NVECTOR) -DNDIM=2 -DNPRE=$(NPRE) -DNVAR=$(NVAR) -DSOLVER$(SOLVER)
DEFINES3D = -DNVECTOR=$(NVECTOR) -DNDIM=3 -DNPRE=$(NPRE) -DNVAR=$(NVAR) -DSOLVER$(SOLVER)
FFLAGS = $(FCFLAGS) -O3 -frecord-marker=4 -fbacktrace -ffree-line-length-none -g 
#FFLAGS = -O3 -frecord-marker=4 -fbacktrace -ffree-line-length-none -g 
FFLAGS += -x f95-cpp-input -DWITHOUTMPI
FFLAGS1D = $(FFLAGS) $(DEFINES1D)
FFLAGS2D = $(FFLAGS) $(DEFINES2D)
FFLAGS3D = $(FFLAGS) $(DEFINES3D)


BUILDDIR = build
BUILDDIR2D = build_2d
BUILDDIR1D = build_1d

OBJS = $(BUILDDIR)/interface.o
OBJS1D = $(BUILDDIR1D)/interface.o
OBJS2D = $(BUILDDIR2D)/interface.o

LIBNAME = "libramses.a"
CODELIB = $(BUILDDIR)/bin/$(LIBNAME)
CODELIB1D = $(BUILDDIR1D)/bin/$(LIBNAME)
CODELIB2D = $(BUILDDIR2D)/bin/$(LIBNAME)

CLASSNAME=RamsesInterface

all: ramses_worker

$(BUILDDIR)/bin/Makefile: $(SRC)/Makefile
	cp -r $(SRC)/.. $(BUILDDIR)
$(BUILDDIR2D)/bin/Makefile: $(SRC)/Makefile
	cp -r $(SRC)/.. $(BUILDDIR2D)
$(BUILDDIR1D)/bin/Makefile: $(SRC)/Makefile
	cp -r $(SRC)/.. $(BUILDDIR1D)

$(CODELIB): $(BUILDDIR)/bin/Makefile
	make -C $(BUILDDIR)/bin all F90="$(FC)" FFLAGS="$(FFLAGS3D)"
$(CODELIB2D): $(BUILDDIR2D)/bin/Makefile
	make -C $(BUILDDIR2D)/bin all F90="$(FC)" FFLAGS="$(FFLAGS2D)"
$(CODELIB1D): $(BUILDDIR1D)/bin/Makefile
	make -C $(BUILDDIR1D)/bin all F90="$(FC)" FFLAGS="$(FFLAGS1D)"

ramses_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) -I$(SRC) $(SC_FLAGS) $(FS_FLAGS) worker_code.f90 -o $@ $(OBJS) $(CODELIB) $(SC_CLIBS) $(FS_LIBS)
ramses_worker_2d: worker_code.f90 $(CODELIB2D) $(OBJS2D)
	$(MPIFC) -I$(SRC) $(SC_FLAGS) $(FS_FLAGS) worker_code.f90 -o $@ $(OBJS2D) $(CODELIB2D) $(SC_CLIBS) $(FS_LIBS)
ramses_worker_1d: worker_code.f90 $(CODELIB1D) $(OBJS1D)
	$(MPIFC) -I$(SRC) $(SC_FLAGS) $(FS_FLAGS) worker_code.f90 -o $@ $(OBJS1D) $(CODELIB1D) $(SC_CLIBS) $(FS_LIBS)

$(OBJS): interface.f90 $(CODELIB)
	$(FC) -c -I$(BUILDDIR)/bin $(FFLAGS3D) $< -o $@
$(OBJS2D): interface.f90 $(CODELIB2D)
	$(FC) -c -I$(BUILDDIR2D)/bin $(FFLAGS2D) $< -o $@
$(OBJS1D): interface.f90 $(CODELIB1D)
	$(FC) -c -I$(BUILDDIR1D)/bin $(FFLAGS1D) $< -o $@

worker_code.f90: interface.py 
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

clean:
	rm -f *.o *.pyc
	rm -f worker_code.f90
	rm -f ramses_worker ramses_worker_mpi
	rm -Rf $(BUILDDIR)
	rm -Rf $(BUILDDIR2D)
	rm -Rf $(BUILDDIR1D)
	make -C $(SRC) clean
