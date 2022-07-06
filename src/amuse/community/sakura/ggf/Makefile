# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk

PYTHON ?= python

all: code sakura_worker

code: install

ifdef DOWNLOAD_CODES
install:
	make -C . install_tupan
else
install:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set."
	@echo "Sakura needs Tupan to be installed."
	@echo "Tupan will not be downloaded and build."
	@echo "If you do want Tupan, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif 

install_tupan:
	(cd ../tupan; make DOWNLOAD_CODES=1)


sakura_worker: interface.py
	$(CODE_GENERATOR) --type=py --mode=mpi -x amuse.community.sakura.interface SakuraInterface SakuraImplementation -o $@
	

clean:
	rm -f *.pyc;
	rm -f sakura_worker sakura_worker_sockets

distclean: clean
	(cd ../tupan; make distclean)

