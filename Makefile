-include config.mk

PYTHON ?= python2.6
VERSION ?= undefined
CLEAN ?= yes

CONFIGURE_ERROR=

export PYTHONPATH := $(PYTHONPATH):$(PWD)/src:$(PWD)/test

all: build.py
	@-mkdir -p test_results
	$(PYTHON) setup.py generate_main
	$(PYTHON) setup.py build_codes --inplace

build.py:
	$(error the code is not configured, please run configure first)

allinbuild:
	$(PYTHON) setup.py build

docclean:
	make -C doc clean

clean:
	$(PYTHON) setup.py clean
	$(PYTHON) setup.py clean_codes --inplace

distclean:
	-rm -f config.mk
	-rm -f support/config.py
	-rm -f support/config.pyc
	-rm -f src/amuse/config.py
	-rm -f src/amuse/config.pyc
	-rm -f amuse.sh
	-rm -f iamuse.sh
	-rm -f ibis-deploy.sh
	-rm -f build.py
	
	-rm -f test/*.000 test/fort.* test/perr test/pout test/test.h5 test/*.log
	-rm -f test/codes_tests/perr test/codes_tests/pout
	-rm -f test/core_tests/plummer_back_100.ini
	-rm -f test/test_python_implementation test/twobody
	
	$(PYTHON) setup.py clean
	$(PYTHON) setup.py dist_clean
	$(PYTHON) setup.py dist_clean --inplace
	
	make -C doc clean
	-find src -name "*.pyc" -exec rm \{} \;
	-find src -type d -name "__pycache__" -exec rm -Rf \{} \;
	-find src -type d -name "ccache" -exec rm -Rf \{} \;
	-rm -Rf build

tests:
	$(PYTHON) setup.py tests
	

doc:
	$(PYTHON) setup.py -q build_latex

html:
	make -C doc html

latexpdf:
	make -C doc latexpdf

ctags:
	find src -name "*.py" | xargs ctags
	find src -name "*.cc" | xargs ctags -a
	find src -name "*.[cCfFhH]" | xargs ctags -a
	find src -name "*.cpp" | xargs ctags -a

release: distclean
	sed 's/version = .*/version = "$(VERSION)",/' < setup.py > releasesetup.py
	cp setup.py setup.py.bck
	mv releasesetup.py setup.py 
	make -C doc release
	python setup.py  -q sdist
	cp setup.py.bck setup.py

nightly:
	make -C doc release
	sed 's/version = .*/version = "$(VERSION)",/' < setup.py > nightlysetup.py
	cp setup.py setup.py.bck
	mv nightlysetup.py setup.py 
	python setup.py sdist
	cp setup.py.bck setup.py
	
debian:
	$(PYTHON) ./support/debian.py


%.code:
ifneq (,$(findstring s,$(MAKEFLAGS)))
	$(PYTHON) setup.py build_code --inplace --clean=$(CLEAN) --code-name=$*
else
	$(PYTHON) setup.py -v build_code --inplace --clean=$(CLEAN) --code-name=$*
endif

