
PYTHON=python

all:
	$(PYTHON) setup.py build

tests:
	$(PYTHON) setup.py tests

doc:
	$(PYTHON) setup.py build_latex
