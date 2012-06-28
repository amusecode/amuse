AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC  ?= mpif90
MPIF90 ?= $(MPIFC)
PYTHON ?= python
FC      = $(MPIFC)

FCLAGS   += -fno-automatic
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libmmc.a

MMC_DATA_DIR = ../../../../data/mmc
CODE_GENERATOR = $(AMUSE_DIR)/build.py

SRCINC = -I./src

all: mmc_worker

clean:
	$(RM) -f *.o *.pyc worker_code.f90  
	$(RM) -f *~ worker_code worker_code.f90 interface.f90
	$(RM) -f *.tmp
	$(RM) $(MMC_DATA_DIR)/*.dat
	$(RM) -f mmc.mod
	$(RM) -f src/amuse_interface2.f
	make -C src clean
	make -C amuselib clean

distclean: clean
	$(RM) -f mmc.mod
	
$(CODELIB): src/amuse_interface2.f
	make -C src
	make -C amuselib

interface.f90 src/amuse_interface2.f: makeinterface.py
	$(PYTHON) makeinterface.py

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py mmcInterface -o $@

mmc_worker: worker_code.f90 $(OBJS) $(CODELIB) 
	$(MPIFC) $(FCFLAGS) $(OBJS) $(SRCINC) $< $(CODELIB) -o $@

interface.o: interface.f90 
	$(FC) $(FCLAGS) $(SRCINC) -c -o $@ $<

