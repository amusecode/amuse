ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif

-include ${AMUSE_DIR}/config.mk 				  # NG has some doubts

MPIF90 ?= $(MPIFC) 	
	 							 
FORTRAN ?= gfortran						   

F77FLAGS =  									  # 
ifneq (,$(findstring xlf, $(notdir $(FORTRAN))))  # 
F77FLAGS += -qfixed  							  # 
endif  											  # NG has some doubts


CODE_GENERATOR ?= $(AMUSE_DIR)/build.py

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

mosse_worker:	worker_code.f90  interface.o  $(OBJ)
	make -C src mosse CMPLR=$(FORTRAN) FFLAGS="$(FCFLAGS) $(F77FLAGS)"
	$(MPIF90) $(FS_FLAGS) $(LDFLAGS) $^ -o $@ $(FS_LIBS) $(LIBS)
	
worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MOSSEInterface -o $@
	
.f.o: $<                                                                 
	$(FORTRAN) -c $(F77FLAGS) $(FCFLAGS) -o $@ $< 
	
%.o: %.f90
	$(FORTRAN) $(FFLAGS)  $(FCFLAGS) -c -o $@ $<

