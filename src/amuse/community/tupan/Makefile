ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
tupan.tar.gz:
	$(DOWNLOAD) https://github.com/ggf84/tupan/archive/refs/heads/master.tar.gz >$@

src/tupan: tupan.tar.gz
	tar xf $<
	mv tupan-master src/tupan
	# TODO: does the configured directory work?
	cp tupan.cfg src/tupan/


# Building and installing the tupan Python package

# Ideally, Tupan would have Conda and PyPI packages, so that we could just have it as a
# top-level dependency, and it would also help if it had its dependencies in the
# setup.py metadata, but this works. Note that we don't install everything in Tupan's
# requirements.txt, because not all of them are needed for use with AMUSE.

.PHONY: install_tupan
install_tupan: | src/tupan
ifneq ($(CONDA_PREFIX),)
	conda install --yes cffi numpy scipy matplotlib
else ifneq ($(VIRTUAL_ENV),)
	pip install cffi numpy scipy matplotlib
endif
	cd src/tupan && pip install .


# Building the workers
tupan_worker: interface.py
	PYTHONPATH=. amusifier --type=py --mode=mpi -x interface TupanInterface TupanImplementation -o $@


# Which packages contain which workers?
amuse-tupan_contains: tupan_worker


# Building and installing packages
develop-%: %_contains install_tupan
	support/shared/uninstall.sh $*
	python -m pip install -e packages/$*

install-%: %_contains install_tupan
	support/shared/uninstall.sh $*
	python -m pip install packages/$*

package-%: %_contains
	python3 -m pip install -vv --no-cache-dir --no-deps --no-build-isolation --prefix ${PREFIX} packages/$*

test-%:
	cd packages/$* && pytest


# Cleaning up
.PHONY: clean
clean:
	rm -rf *worker*

.PHONY: distclean
distclean: clean
	rm -rf tupan.tar.gz src/tupan
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
