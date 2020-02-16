ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk

PYTHON ?= python

ARGPARSE_AVAILABLE := $(shell $(PYTHON) -c "import argparse"  1>&2 2> /dev/null && echo "yes" || echo "no")

CODE_GENERATOR ?= $(AMUSE_DIR)/build.py

CLASSNAME=KromeInterface

FFLAGS += $(FCFLAGS) -O2 -ffree-line-length-none
INCL = -I./src/ -I./src/krome/build/

all: src/krome src/libchem.a test krome_worker

ifdef DOWNLOAD_CODES
src/krome:
	make -C . download
else
src/krome:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. Krome will not be downloaded and build."
	@echo "If you do want Athena, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif 

latest_tested.tar.gz: URL=$(shell wget "http://kromepackage.org/test/get.php?ref=tested_url" -q -O -)
latest_tested.tar.gz:
	wget $(URL) -O latest_tested.tar.gz

archive_root: ARCHIVE_ROOT=$(shell tar -tf latest_tested.tar.gz | head -1 | sed -e 's/\/.*//')
archive_root:
	tar xvzf latest_tested.tar.gz
	rm -rf src/krome
	mv $(ARCHIVE_ROOT) src/krome

download: latest_tested.tar.gz  archive_root

update: distclean download

download_git:
	git clone https://bitbucket.org/tgrassi/krome.git
	mv krome src/krome

test:
	@echo
	@echo "Testing import of modules required for krome (argparse):"
ifeq ($(ARGPARSE_AVAILABLE),no)
	$(error "Python imports not available: argparse")
endif
	@echo "Tests successful!"
	@echo
	make -C src/ test

src/libchem.a: src/krome
	make -C src/ libchem.a FFLAGS="$(FFLAGS)" FC="$(FC)"

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 $< $(CLASSNAME) -o $@

krome_worker: src/libchem.a worker_code.f90 interface.o
	$(MPIFC) $(FFLAGS) $(INCL) $(SC_FLAGS) $(FS_FLAGS) $(LDFLAGS) worker_code.f90 interface.o -o $@ src/libchem.a  $(LIBS) $(SC_FCLIBS) $(FS_LIBS)  $(LIBS)


%.o: %.f90
	$(FC) $(FFLAGS) $(INCL) -c -o $@ $<

clean:
	make -C src/ clean
	rm -f *.pyc
	rm -f interface.o krome_worker.f90 worker_code.f90
	rm -f krome_worker 

distclean: clean
	make -C src/ distclean 
	rm -rf src/krome
	rm -rf latest_tested.tar.gz

