# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx

OBJ = interface.o

CODE_GENERATOR = $(PYTHON) $(AMUSE_DIR)/build.py

all: fastkick_worker

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py FastKickInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py FastKickInterface -o $@

fastkick_worker: worker_code.cc worker_code.h $(OBJ)
	$(MPICXX) $(CXXFLAGS) $< $(OBJ) -o $@ 

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 

clean:
	rm -f *.o *.pyc worker_code.cc worker_code.h
	rm -f fastkick_worker
