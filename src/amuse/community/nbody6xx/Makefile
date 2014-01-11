# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
FC      = $(MPIFC)

LDFLAGS  += -lm $(MUSE_LD_FLAGS) 

OBJS = interface.o

CODELIB = src/libnbody6xx.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

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
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $< $(OBJS) $(CODELIB) $(FS_LIBS) -o $@  $(LIBS)

%.o: %.f90 
	$(FC) $(FCFLAGS) -Isrc/ -ffixed-form -xf77-cpp-input -c -o $@ $<
