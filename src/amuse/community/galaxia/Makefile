# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(strip $(AMUSE_DIR)),)
    AMUSE_DIR = ../../../..
endif
-include ${AMUSE_DIR}/config.mk

ifeq ($(strip $(MPIFC)),)
    MPIFC = mpif90
endif
ifeq ($(strip $(FC)),)
    FC = gfortran
endif

ifeq ($(strip $(FSLIB)),)
    FSLIB = -L$(AMUSE_DIR)/lib/forsockets -lforsockets 
endif
ifeq ($(strip $(FSINC)),)
    FSINC = -I$(AMUSE_DIR)/lib/forsockets
endif

CODE_GENERATOR = $(AMUSE_DIR)/build.py

OBJ	= interface.o

all:    GalaxyModel_worker

clean:
	rm -f  *~ *.pyc *.mod *.o GalaxyModel_worker worker_code.f90


$(OBJ): interface.f90
	$(FC) -c $(FCFLAGS) -llapack -o $@ $< 

GalaxyModel_worker:	worker_code.f90  $(OBJ)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $^ -o $@  $(FS_LIBS) -llapack

.f90.o: $<
	$(FC) -c $(FCFLAGS)  -llapack -o $@ $< 

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py  BarAndSpiralsInterface -o $@

