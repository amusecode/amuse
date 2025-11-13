ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk


CLASSNAME=UCLchemInterface

FFLAGS += $(FCFLAGS) -std=legacy -O2 -ffree-line-length-none
INCL = -I./src/ 
export FFLAGS

all: uclchem_worker

.PHONY: src/libchem.a
src/libchem.a: 
	$(MAKE) -C src/ libchem.a 

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

uclchem_worker: src/libchem.a worker_code.f90 interface.o
	$(MPIFC) $(FFLAGS) $(INCL) $(FS_FLAGS) $(LDFLAGS) worker_code.f90 interface.o -o $@ src/libchem.a  $(LIBS) $(FS_LIBS)  $(LIBS)

interface.o: src/libchem.a

%.o: %.f90
	$(MPIFC) $(FFLAGS) $(INCL) -c -o $@ $<

clean:
	make -C src/ clean
	rm -f *.pyc
	rm -f interface.o uclchem_worker.f90 worker_code.f90
	rm -f uclchem_worker 

distclean: clean

	

