
PYTHON=python

export PYTHONPATH := $(PWD)/src:$(PWD)/test	

all:
	$(PYTHON) setup.py build

docclean:
	make -C doc clean
	
clean: docclean
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

%.code:
	$(PYTHON) setup.py code --code-name=$*
