ifndef AMUSE_DIR
	AMUSE_DIR=../../../..
endif

-include ${AMUSE_DIR}/config.mk 				  # NG has some doubts

MPIFC ?= mpif90  								  # 
FORTRAN ?= gfortran						  		  # 

F77FLAGS =  									  # 
ifneq (,$(findstring xlf, $(notdir $(FORTRAN))))  # 
F77FLAGS += -qfixed  							  # 
endif  											  # NG has some doubts


CODE_GENERATOR = $(AMUSE_DIR)/build.py

OBJ	= src/comenv.o src/corerd.o src/deltat.o src/dgcore.o \
	  src/evolv3.o src/gntage.o src/hrdiag.o src/instar.o \
	  src/kick.o src/mix.o src/mlwind.o src/mrenv.o \
	  src/ran3.o src/rl.o \
	  src/pisn.f src/eddington.f src/fallback.f \
	  src/star.o src/zcnsts.o src/zfuncs.o
	  
MOBSEOBJ = src/comenv.o src/corerd.o src/deltat.o src/dgcore.o \
         src/evolv3.o src/gntage.o src/hrdiag.o src/instar.o \
         src/kick.o src/mix.o src/mlwind.o src/mrenv.o \
         src/ran3.o src/rl.o \
		 src/pisn.f src/eddington.f src/fallback.f \
		 src/star.o src/zcnsts.o src/zfuncs.o

all:	mobse_worker

clean:
	$(RM) -f  *~ *.so *.o src/*.o 
	$(RM) -f src/mobse mobse_worker worker_code worker_code.f90 
	$(RM) -f worker_code-sockets.f90 mobse_worker_sockets
	$(RM) -f amuse_worker_module.mod
	$(RM) -f forsockets.mod	

mobse_worker:	worker_code.f90  interface.o  $(MOBSEOBJ)         
	make -C src mobse CMPLR=$(FORTRAN) FFLAGS="$(FFLAGS) $(F77FLAGS)"
	$(MPIFC) $(FFLAGS) $(FS_FLAGS) $(LDFLAGS) $^ -o $@ $(FS_LIBS)  $(LIBS)

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MOBSEInterface -o $@

.f.o: $<
	$(FORTRAN) -c $(F77FLAGS) $(FFLAGS) -o $@ $< 
