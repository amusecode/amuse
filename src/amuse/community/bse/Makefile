
MPIF90 ?= mpif90
FORTRAN ?= gfortran
AMUSE_DIR ?=../../../..

FSOPTS ?= -L$(AMUSE_DIR)/lib/forsockets -lforsockets -I$(AMUSE_DIR)/lib/forsockets

CODE_GENERATOR = $(AMUSE_DIR)/build.py

OBJ	= src/comenv.o src/corerd.o src/deltat.o src/dgcore.o \
	  src/evolv1.o src/evolv2.o src/gntage.o src/hrdiag.o \
	  src/instar.o src/kick.o src/mix.o src/mlwind.o \
	  src/mrenv.o src/popbin.o src/ran3.o src/rl.o \
	  src/sse.o src/star.o src/zcnsts.o src/zfuncs.o
	  
BSEOBJ = src/comenv.o src/corerd.o src/deltat.o src/dgcore.o \
         src/evolv2.o src/gntage.o src/hrdiag.o src/instar.o \
         src/kick.o src/mix.o src/mlwind.o src/mrenv.o \
         src/ran3.o src/rl.o src/star.o src/zcnsts.o src/zfuncs.o

all:    bse_worker

clean:
	/bin/rm -f *~ *.so *.o src/*.o src/bse bse_worker worker_code worker_code.f90 worker_code-sockets.f90 bse_worker_sockets

bse_worker:	worker_code.f90  BSE_code.o  $(BSEOBJ)
	make -C src bse CMPLR=$(FORTRAN) FFLAGS="$(FFLAGS)"
	$(MPIF90) $(FFLAGS) $^ -o $@ 
	
bse_worker_sockets:	worker_code-sockets.f90  BSE_code.o  $(BSEOBJ)
	make -C src bse CMPLR=$(FORTRAN) FFLAGS="$(FFLAGS)"
	$(MPIF90) $(FFLAGS) $(FSOPTS) $^ -o $@ 

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py BSEInterface -o $@
	
worker_code-sockets.f90: interface.py
	$(CODE_GENERATOR) --type=f90 --mode=sockets interface.py BSEInterface -o $@

.f.o: $<
	$(FORTRAN) -c $(FFLAGS) -o $@ $< 
