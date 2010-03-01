
FORTRAN ?= gfortran
AMUSE_DIR ?= ../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

OBJ	= src/deltat.o src/evolv1.o src/hrdiag.o src/kick.o \
	  src/mlwind.o src/mrenv.o src/ran3.o src/star.o \
	  src/zcnsts.o src/zfuncs.o

all:    worker_code


clean:
	/bin/rm -f *~ *.so *.o src/*.o src/sse worker_code worker_code.f90

worker_code:	worker_code.f90  SSE_muse_interface.o  $(OBJ)
	make -C src sse CMPLR=$(FORTRAN)
	mpif90 $^ -o $@ 


worker_code.f90: muse_stellar_mpi.py
	$(CODE_GENERATOR) --type=f90 muse_stellar_mpi.py SSE -o $@

.f.o: $<                                                                 
	$(FORTRAN) -c $(FFLAGS) -o $@ $< 
