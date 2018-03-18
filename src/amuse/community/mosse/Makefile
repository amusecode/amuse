ifndef AMUSE_DIR
	AMUSE_DIR=../../../..
endif

-include ${AMUSE_DIR}/config.mk 				  # NG has some doubts

MPIF90 ?= $(MPIFC) 	
	 							 
FORTRAN ?= gfortran						   

F77FLAGS =  									  # 
ifneq (,$(findstring xlf, $(notdir $(FORTRAN))))  # 
F77FLAGS += -qfixed  							  # 
endif  											  # NG has some doubts


CODE_GENERATOR = $(AMUSE_DIR)/build.py

FSLIBS ?= -L$(AMUSE_DIR)/lib/forsockets -lforsockets
FSINC ?= -I$(AMUSE_DIR)/lib/forsockets

AM_LIBS = -L$(AMUSE_DIR)/lib/amuse_mpi -lamuse_mpi
AM_CFLAGS = -I$(AMUSE_DIR)/lib/amuse_mpi

OBJ	= src/deltat.o src/evolv1.o src/hrdiag.o  \
	  src/kick.o src/mlwind.o src/mrenv.o \
	  src/ran3.o src/pisn.f src/eddington.f src/fallback.f \
	  src/star.o src/zcnsts.o src/zfuncs.o
	  
all:	mosse_worker

clean:
	$(RM) -f *~ *.so *.pyc *.o src/*.o src/mosse worker_code worker_code.f90
	$(RM) -f mosse_worker mosse_worker_sockets worker_code-sockets.f90
	$(RM) -f amuse_worker_module.mod
	$(RM) -f forsockets.mod
	$(RM) -f mosse_cython.* mosse_worker_* mosse_c.* cinterfacemodule.mod

mosse_worker:	worker_code.f90  interface.o  $(OBJ)
	make -C src mosse CMPLR=$(FORTRAN) FFLAGS="$(FCFLAGS) $(F77FLAGS)"
	$(MPIF90) $(FS_FLAGS) $(LDFLAGS) $^ -o $@ $(FS_LIBS) $(LIBS)
	
mosse_worker_cython: mosse_cython.so
	$(CODE_GENERATOR) --type=cython -m script -x amuse.community.mosse.interface MOSSEInterface -o $@ --cython-import mosse_cython
	
mosse_cython.so: mosse_cython.o $(CODELIB) $(OBJ) interface.o mosse_c.o
	$(MPIFC) -shared $(CXXFLAGS) $(PYTHONDEV_LDFLAGS) $(AM_CFLAGS) $(SC_FLAGS) $(LDFLAGS) $(OBJ) mosse_cython.o interface.o  -o $@ $(SC_CLIBS) $(AM_LIBS) mosse_c.o

mosse_cython.o: mosse_cython.c
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(AM_CFLAGS) $(PYTHONDEV_CFLAGS) -c -o $@ $< 

mosse_cython.c: mosse_cython.pyx 
	$(CYTHON) $< -o $@

mosse_c.f90: interface.py
	$(CODE_GENERATOR) --type=cython -m interface -x amuse.community.mosse.interface MOSSEInterface -o $@ --prefix aaa_

mosse_cython.pyx: interface.py
	$(CODE_GENERATOR) --type=cython -m mpi interface.py MOSSEInterface -o $@ --prefix aaa_


worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MOSSEInterface -o $@
	
.f.o: $<                                                                 
	$(FORTRAN) -c $(F77FLAGS) $(FCFLAGS) -o $@ $< 
	
%.o: %.f90
	$(FORTRAN) $(FFLAGS)  $(FCFLAGS) -c -o $@ $<

