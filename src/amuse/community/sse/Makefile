# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPIF90 ?= $(MPIFC)

FORTRAN ?= gfortran

F77FLAGS = 
ifneq (,$(findstring xlf, $(notdir $(FORTRAN))))
F77FLAGS += -qfixed
endif

CODE_GENERATOR = $(AMUSE_DIR)/build.py

FSLIBS ?= -L$(AMUSE_DIR)/lib/forsockets -lforsockets
FSINC ?= -I$(AMUSE_DIR)/lib/forsockets


AM_LIBS = -L$(AMUSE_DIR)/lib/amuse_mpi -lamuse_mpi
AM_CFLAGS = -I$(AMUSE_DIR)/lib/amuse_mpi

OBJ	= src/deltat.o src/evolv1.o src/hrdiag.o src/kick.o \
	  src/mlwind.o src/mrenv.o src/ran3.o src/star.o \
	  src/zcnsts.o src/zfuncs.o

all:    sse_worker


clean:
	$(RM) -f *~ *.so *.pyc *.o src/*.o src/sse worker_code worker_code.f90
	$(RM) -f sse_worker sse_worker_sockets worker_code-sockets.f90
	$(RM) -f amuse_worker_module.mod
	$(RM) -f forsockets.mod
	$(RM) -f sse_cython.* sse_worker_* sse_c.* cinterfacemodule.mod

sse_worker:	worker_code.f90  interface.o  $(OBJ)
	make -C src sse CMPLR=$(FORTRAN) FFLAGS="$(FCFLAGS) $(F77FLAGS)"
	$(MPIF90) $(FS_FLAGS) $(LDFLAGS) $^ -o $@ $(FS_LIBS) $(LIBS)
	
sse_worker_cython: sse_cython.so
	$(CODE_GENERATOR) --type=cython -m script -x amuse.community.sse.interface SSEInterface -o $@ --cython-import sse_cython
	
sse_cython.so: sse_cython.o $(CODELIB) $(OBJ) interface.o sse_c.o
	$(MPIFC) -shared $(CXXFLAGS) $(PYTHONDEV_LDFLAGS) $(AM_CFLAGS) $(SC_FLAGS) $(LDFLAGS) $(OBJ) sse_cython.o interface.o  -o $@ $(SC_CLIBS) $(AM_LIBS) sse_c.o

sse_cython.o: sse_cython.c
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(AM_CFLAGS) $(PYTHONDEV_CFLAGS) -c -o $@ $< 

sse_cython.c: sse_cython.pyx 
	$(CYTHON) $< -o $@

sse_c.f90: interface.py
	$(CODE_GENERATOR) --type=cython -m interface -x amuse.community.sse.interface SSEInterface -o $@ --prefix aaa_

sse_cython.pyx: interface.py
	$(CODE_GENERATOR) --type=cython -m mpi interface.py SSEInterface -o $@ --prefix aaa_


worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py SSEInterface -o $@
	
.f.o: $<                                                                 
	$(FORTRAN) -c $(F77FLAGS) $(FCFLAGS) -o $@ $< 
	
%.o: %.f90
	$(FORTRAN) $(FFLAGS)  $(FCFLAGS) -c -o $@ $<

distclean: clean
