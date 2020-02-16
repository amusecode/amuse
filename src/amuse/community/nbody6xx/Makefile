# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
FORTRAN ?= $(FC)

LDFLAGS  += -lm $(MUSE_LD_FLAGS) 

ifeq ($(findstring gfortran, $(notdir $(FORTRAN))), gfortran)
FCFLAGS2+=-ffixed-form -xf77-cpp-input
endif

ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
FCFLAGS2=-fpp -fixed
endif

OBJS = interface.o

CODELIB = src/libnbody6xx.a

CODE_GENERATOR ?= $(AMUSE_DIR)/build.py

all: nbody6xx_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) -f *~ nbody6xx_worker worker_code.f90
	$(RM) -f amuse_interface.mod 
	make -C src clean

$(CODELIB):
ifeq ($(MPI_ENABLED), no)
	make -C src amuse
else
	make -C src amusepar
endif

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py Nbody6xxInterface -o $@

nbody6xx_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $(LDFLAGS) $< $(OBJS) $(CODELIB) $(FS_LIBS) -o $@  $(LIBS)

%.o: %.f90 
	$(FC) $(FCFLAGS) $(FCFLAGS2) -Isrc/ -c -o $@ $<
