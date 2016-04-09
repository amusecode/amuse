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

CODE_GENERATOR = $(AMUSE_DIR)/build.py

AM_LIBS = -L$(AMUSE_DIR)/lib/amuse_mpi -lamuse_mpi
AM_CFLAGS = -I$(AMUSE_DIR)/lib/amuse_mpi


all: aarsethzare_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ aarsethzare_worker worker_code.f90 aarsethzare_worker_sockets worker_code-sockets.f90
	$(RM) -f *.mod
	$(RM) -f aarsethzare_worker_cython aarsethzare_cython.* aarsethzare_c.*
	make -C src clean

$(CODELIB):  src/*.f
	make -C src all

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py AarsethZareInterface -o $@
	
aarsethzare_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $(LDFLAGS) -cpp worker_code.f90 $(OBJS) $(CODELIB) $(FS_LIBS) -o $@  $(LIBS)
	
aarsethzare_worker_cython: aarsethzare_cython.so
	$(CODE_GENERATOR) --type=cython -m script -x amuse.community.aarsethzare.interface AarsethZareInterface -o $@ --cython-import aarsethzare_cython
	
aarsethzare_cython.so: aarsethzare_cython.o $(CODELIB) $(OBJS) aarsethzare_c.o
	$(MPIFC) -shared $(CXXFLAGS) $(PYTHONDEV_LDFLAGS) $(AM_CFLAGS) $(SC_FLAGS) $(LDFLAGS) aarsethzare_cython.o  -o $@ $(SC_CLIBS) $(AM_LIBS) $(OBJS) $(CODELIB) aarsethzare_c.o

aarsethzare_cython.o: aarsethzare_cython.c
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(AM_CFLAGS) $(PYTHONDEV_CFLAGS) -c -o $@ $< 

aarsethzare_cython.c: aarsethzare_cython.pyx 
	$(CYTHON) $< -o $@

aarsethzare_c.f90: interface.py
	$(CODE_GENERATOR) --type=cython -m interface -x amuse.community.aarsethzare.interface AarsethZareInterface -o $@ --prefix aaa_

aarsethzare_cython.pyx: interface.py
	$(CODE_GENERATOR) --type=cython -m mpi interface.py AarsethZareInterface -o $@ --prefix aaa_

interfacec.f90: interface.py
	$(CODE_GENERATOR) --type=cython -m interface -x amuse.community.aarsethzare.interface AarsethZareInterface -o $@
	
worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py AarsethZareInterface -o $@

	
%.o: %.f90
	$(FORTRAN) $(FFLAGS)  $(FCFLAGS) -c -o $@ $<
