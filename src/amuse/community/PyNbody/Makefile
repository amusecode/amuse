all: build install

build:
	(cd PyNbody; python setup.py build)

install:
	(cd PyNbody; python setup.py install --user)

clean:
	rm -f *.pyc; (cd PyNbody; python setup.py clean -a)

