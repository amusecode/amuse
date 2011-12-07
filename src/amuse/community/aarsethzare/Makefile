# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPIFC ?= mpif90
FORTRAN ?= gfortran
FC ?= $(FORTRAN)

FFLAGS   += -g
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libaarsethzare.a

FSOPTS ?= -L$(AMUSE_DIR)/lib/forsockets -lforsockets -I$(AMUSE_DIR)/lib/forsockets

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: aarsethzare_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ aarsethzare_worker worker_code.f90 aarsethzare_worker_sockets worker_code-sockets.f90
	$(RM) -f forsockets.mod
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py AarsethZareInterface -o $@
	
worker_code-sockets.f90: interface.py
	$(CODE_GENERATOR) --type=f90 --mode=sockets interface.py AarsethZareInterface -o $@

aarsethzare_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) worker_code.f90 $(OBJS) $(CODELIB) -o $@
	
aarsethzare_worker_sockets: worker_code-sockets.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(FSOPTS) worker_code-sockets.f90 $(OBJS) $(CODELIB) -o $@ 
	

%.o: %.f90
	$(FORTRAN) $(FFLAGS) -c -o $@ $<
