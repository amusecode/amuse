ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

MPIFC  ?= mpif90
MPIF90 ?= $(MPIFC)
PYTHON ?= python
FC      = $(MPIFC)

FCFLAGS   += -fno-automatic
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libmmc.a

MMC_DATA_DIR = ../../../../data/mmc

SRCINC = -I./src

all: mmc_worker

clean:
	$(RM) -f *.o *.pyc worker_code.f90
	$(RM) -f *~ worker_code worker_code.f90 interface.f90
	$(RM) -f *.tmp mmc_worker
	$(RM) -f $(MMC_DATA_DIR)/*.dat
	$(RM) -f mmc.mod
	$(RM) -f src/amuse_interface2.f
	-$(MAKE) -C src clean
	-$(MAKE) -C amuselib clean

distclean: clean
	$(RM) -f mmc.mod
	
$(CODELIB): src/amuse_interface2.f
	$(MAKE) -C src
	$(MAKE) -C amuselib

interface.f90 src/amuse_interface2.f: makeinterface.py
	$(PYTHON) makeinterface.py

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py mmcInterface -o $@

mmc_worker: worker_code.f90 $(OBJS) $(CODELIB)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $(LDFLAGS) $(OBJS) $(SRCINC) $< $(CODELIB) $(FS_LIBS) -o $@  $(LIBS)

interface.o: interface.f90
	$(FC) $(FCFLAGS) $(SRCINC) -c -o $@ $<

