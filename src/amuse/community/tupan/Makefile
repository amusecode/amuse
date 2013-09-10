AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

PYTHON ?= python

PYOPENGL_AVAILABLE := $(shell $(PYTHON) -c "import OpenGL"  1>&2 2> /dev/null && echo "yes" || echo "no")
PYIMAGE_AVAILABLE := $(shell $(PYTHON) -c "import Image"  1>&2 2> /dev/null && echo "yes" || echo "no")
PYOPENCL_AVAILABLE := $(shell $(PYTHON) -c "import pyopencl"  1>&2 2> /dev/null && echo "yes" || echo "no")

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: code build install test tupan_worker

code: tupan/setup.py

ifdef DOWNLOAD_CODES
tupan/setup.py:
	make -C . download
else
tupan/setup.py:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. tupan will not be downloaded and build."
	@echo "If you do want Athena, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif 

download:
	git clone -b master https://github.com/ggf84/tupan.git

build:
	(cd tupan; python setup.py build)

install:
	(cd tupan; python setup.py install --user)

test:
	@echo
	@echo "Testing import of modules required for tupan (OpenGL, PIL, pyopencl):"
ifeq ($(PYOPENGL_AVAILABLE),no)
	$(error "Python imports not available: OpenGL")
endif
ifeq ($(PYIMAGE_AVAILABLE),no)
	$(error "Python imports not available: Image")
endif
ifeq ($(PYOPENCL_AVAILABLE),no)
	$(error "Python imports not available: pyopencl")
endif
	@echo "Tests successful!"
	@echo


tupan_worker: interface.py
	$(CODE_GENERATOR) --type=py --mode=mpi -x amuse.community.tupan.interface TupanInterface TupanImplementation -o $@
	

clean:
	rm -f *.pyc; (cd tupan; python setup.py clean -a)
	rm -f tupan_worker tupan_worker_sockets

distclean: clean
	rm -f *.pyc; rm -rf tupan

