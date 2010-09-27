PYTHON=python2.6

export PYTHONPATH := $(PYTHONPATH)/:$(PWD)/src:$(PWD)/test

all:
	@-mkdir -p test_results
	$(PYTHON) setup.py build

docclean:
	make -C doc clean

clean:
	$(PYTHON) setup.py clean

tests:
	$(PYTHON) setup.py tests

doc:
	$(PYTHON) setup.py build_latex

html:
	make -C doc html

latexpdf:
	make -C doc latexpdf

ctags:
	find src -name "*.py" | xargs ctags
	find src -name "*.cc" | xargs ctags -a
	find src -name "*.[cCfFhH]" | xargs ctags -a

release:
	make -C doc release
	python setup.py sdist

debian:
	$(PYTHON) ./support/debian.py

%.code:
	$(PYTHON) setup.py code --code-name=$*
