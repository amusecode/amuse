all: download build install

download:
	git clone git://github.com/GuilhermeFerrari/PyNbody.git

build:
	(cd PyNbody; python setup.py build)

install:
	(cd PyNbody; python setup.py install --user)

clean:
	rm -f *.pyc; (cd PyNbody; python setup.py clean -a)

