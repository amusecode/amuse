
PYTHON=python

all:
	$(PYTHON) setup.py build

clean:
	$(PYTHON) setup.py clean


tests:
	$(PYTHON) setup.py tests

doc:
	$(PYTHON) setup.py build_latex
