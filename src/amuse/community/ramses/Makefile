# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPIFC ?= mpif90
FC ?= gfortran

SRC=src/bin

CODE_GENERATOR = $(AMUSE_DIR)/build.py

NOMPILIBS=$(SRC)/libramses.a
#MPILIBS=$(SRC)/libramses_mpi.a

OMPFLAGS= # -openmp -openmp_report0 
PIC=#-fPIC

NVECTOR = 64
NDIM = 1
NPRE = 8
NVAR = 3
SOLVER = hydro
PATCH = 
EXEC = ramses
DEFINES = -DNVECTOR=$(NVECTOR) -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DNVAR=$(NVAR) -DSOLVER$(SOLVER)

ifeq ($(FC), gfortran)
FFLAGS = -O3 -frecord-marker=4 -fbacktrace -ffree-line-length-none -g 
FFLAGS += -x f95-cpp-input $(DEFINES) -DWITHOUTMPI
#FFLAGS += $(DEFINES) -DWITHOUTMPI
endif
#ifeq  ($(FC), ifort)
#FFLAGS += -O -ip -u -g -traceback -fpp1 -DIFORT  $(PIC) $(FCFLAGS)
#endif

MPI_FLAGS = -DMPI # 

FFLAGS+= $(OMPFLAGS) 

CLASSNAME=RamsesInterface

all: ramses_worker

$(NOMPILIBS): $(SRC)/Makefile
#	make -C $(SRC) ramses
	make -C $(SRC) all FFLAGS="$(FFLAGS)" F90="$(FC)"

#$(MPILIBS): $(MPI_BUILDDIR)/Makefile
#	make -C $(MPI_BUILDDIR) amuse_interface_mpi  VPATH=../$(SRC) F90FLAGS1="$(FFLAGS)" FC="$(FC)" MPIFC="$(MPIFC)"
	
ramses_worker: worker_code.f90 interface.o $(NOMPILIBS)
	$(MPIFC) -I$(SRC) $(SC_FLAGS) $(FS_FLAGS) worker_code.f90 -o $@ interface.o $(NOMPILIBS) $(SC_CLIBS) $(FS_LIBS)
#	$(MPIFC) $(FFLAGS) -I$(SRC) $(SC_FLAGS) $(FS_FLAGS) worker_code.f90 -o $@ interface.o $(NOMPILIBS) $(SC_CLIBS) $(FS_LIBS)

#ramses_worker_mpi: $(MPILIBS) worker_code.f90 $(MPI_BUILDDIR)/interface.o
#	$(FC) $(FFLAGS) $(FFLAGS2) $(MPI_FLAGS) -I$(MPI_BUILDDIR) $(SC_FLAGS) $(FS_FLAGS) $^ -o $@ $(MPILIBS) $(SC_MPI_CLIBS) $(FS_LIBS)

interface.o: interface.f90
	$(FC) -c $^ -o $@

worker_code.f90: interface.py 
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

clean:
	rm -f *.o *.pyc
	rm -f worker_code.f90
	rm -f ramses_worker ramses_worker_mpi

veryclean: clean
	make -C $(SRC) clean
