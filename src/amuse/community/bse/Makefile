# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPIFC ?= mpif90
FORTRAN ?= gfortran

F77FLAGS = 
ifneq (,$(findstring xlf, $(notdir $(FORTRAN))))
F77FLAGS += -qfixed
endif


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
	$(RM) -f  *~ *.so *.o src/*.o 
	$(RM) -f src/bse bse_worker worker_code worker_code.f90 
	$(RM) -f worker_code-sockets.f90 bse_worker_sockets
	$(RM) -f amuse_worker_module.mod
	$(RM) -f forsockets.mod
	
bse_worker:	worker_code.f90  interface.o  $(BSEOBJ)
	make -C src bse CMPLR=$(FORTRAN) FFLAGS="$(FFLAGS) $(F77FLAGS)"
	$(MPIFC) $(FFLAGS) $(FS_FLAGS) $^ -o $@ $(FS_LIBS)
	
worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py BSEInterface -o $@
	
.f.o: $<
	$(FORTRAN) -c $(F77FLAGS) $(FFLAGS) -o $@ $< 
