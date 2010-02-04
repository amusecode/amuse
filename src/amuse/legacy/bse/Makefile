
FORTRAN ?= gfortran
AMUSE_DIR ?=../../../..

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

all:    worker_code

clean:
	/bin/rm -f *~ *.so *.o src/*.o src/bse worker_code worker_code.f90

worker_code:	worker_code.f90  BSE_code.o  $(BSEOBJ)
	make -C src bse CMPLR=$(FORTRAN) FFLAGS="$(FFLAGS)"
	mpif90 $(FFLAGS) $^ -o $@ 

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py BSE -o $@

.f.o: $<
	$(FORTRAN) -c $(FFLAGS) -o $@ $< 
