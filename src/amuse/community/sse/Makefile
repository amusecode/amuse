
MPIF90 ?= mpif90

FORTRAN ?= gfortran
AMUSE_DIR ?= ../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

FSLIBS ?= -L$(AMUSE_DIR)/lib/forsockets -lforsockets
FSINC ?= -I$(AMUSE_DIR)/lib/forsockets

OBJ	= src/deltat.o src/evolv1.o src/hrdiag.o src/kick.o \
	  src/mlwind.o src/mrenv.o src/ran3.o src/star.o \
	  src/zcnsts.o src/zfuncs.o

all:    sse_worker


clean:
	rm -f *~ *.so *.o src/*.o src/sse worker_code worker_code.f90
	rm -f sse_worker sse_worker_sockets

sse_worker:	worker_code.f90  SSE_muse_interface.o  $(OBJ)
	make -C src sse CMPLR=$(FORTRAN)
	$(MPIF90) $^ -o $@ 

sse_worker_sockets: worker_code-sockets.f90  SSE_muse_interface.o  $(OBJ)
	make -C src sse CMPLR=$(FORTRAN)
	$(MPIF90) $(FSINC) $^ -o $@ $(FSLIBS) 


worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py SSEInterface -o $@

.f.o: $<                                                                 
	$(FORTRAN) -c $(FFLAGS) -o $@ $< 
