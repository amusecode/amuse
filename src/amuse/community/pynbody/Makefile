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
	python -c "import OpenGL" || (@echo "Python imports not available";$(error Python imports not available: OpenGL))
	python -c "import Image" || (@echo "Python imports not available";$(error Python imports not available: Image))
	python -c "import pyopencl" || (@echo "Python imports not available";$(error Python imports not available: pyopencl))
	@echo "Tests successful!"
	@echo
	
	
clean:
	rm -f *.pyc; (cd PyNbody; python setup.py clean -a)

clean_all:
	rm -f *.pyc; rm -rf PyNbody

