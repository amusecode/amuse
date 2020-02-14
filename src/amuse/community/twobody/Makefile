# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk

CODE_GENERATOR ?= $(AMUSE_DIR)/build.py

all:twobody_worker

twobody_worker: twobody.py
	$(CODE_GENERATOR) --type=py --mode=mpi -x interface TwoBodyInterface TwoBodyImplementation -o $@
	
clean:
	$(RM) -f *.bck *.pyc
	$(RM) -f twobody_worker twobody_worker_sockets
