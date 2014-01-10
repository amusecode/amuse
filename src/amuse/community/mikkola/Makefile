# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPIFC ?= mpif90
FC    = $(MPIFC)

#FFLAGS   += 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

F77FLAGS = 
ifneq (,$(findstring xlf, $(notdir $(FORTRAN))))
F77FLAGS += -qfixed -qextname=flush
endif
export F77FLAGS

CODELIB = src/libMikkola.a

AMUSE_DIR?=../../amuse

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: mikkola_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code worker_code.f90 *.mod mikkola_worker
	make -C src clean

distclean: clean

$(CODELIB):
	make -C src all F77FLAGS="$(F77FLAGS)"

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MikkolaInterface -o $@

mikkola_worker: worker_code.f90 $(CODELIB) $(OBJS)
#	$(MPIFC) $(CXXFLAGS) $(FS_FLAGS) $< $(OBJS) $(CODELIB) -o $@ $(FS_LIBS)
	$(MPIFC) $(FFLAGS) $(FS_FLAGS) $< $(OBJS) $(CODELIB) -o $@ $(FS_LIBS)  $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<
