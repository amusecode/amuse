PYTHON ?= python

PYOPENGL_AVAILABLE := $(shell $(PYTHON) -c "import OpenGL"  1>&2 2> /dev/null && echo "yes" || echo "no")
PYIMAGE_AVAILABLE := $(shell $(PYTHON) -c "import Image"  1>&2 2> /dev/null && echo "yes" || echo "no")
PYOPENCL_AVAILABLE := $(shell $(PYTHON) -c "import pyopencl"  1>&2 2> /dev/null && echo "yes" || echo "no")

all: download build install test

download: PyNbody/setup.py

PyNbody/setup.py:
	git clone git://github.com/GuilhermeFerrari/PyNbody.git

build:
	(cd PyNbody; python setup.py build)

install:
	(cd PyNbody; python setup.py install --user)

test:
	@echo
	@echo "Testing import of modules required for PyNbody (OpenGL, PIL, pyopencl):"
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
	
	
clean:
	rm -f *.pyc; (cd PyNbody; python setup.py clean -a)

clean_all:
	rm -f *.pyc; rm -rf PyNbody

