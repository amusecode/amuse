ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif

# Detecting dependencies
support/config.mk:
	cd support && ./configure

src/tsunami-master: | tsunami-master.tar.gz
	tar xf $<
	mv tsunami-master src/


.PHONY: install_tsunami
install_tsunami: | src/tsunami-master
ifneq ($(CONDA_PREFIX),)
	conda install --yes --channel conda-forge --override-channels numpy "swig=4.2.1" scikit-build-core
endif
	cd src/tsunami-master/ && pip install . -Ccmake.define.spin=ON


tsunami_worker: interface.py install_tsunami
	PYTHONPATH=. amusifier --type=py --mode=mpi -x interface TsunamiInterface TsunamiImplementation -o $@

# Which packages contain which workers?
amuse-tsunami_contains: tsunami_worker


# Building and installing packages
develop-%: %_contains
	support/shared/uninstall.sh $*
	python -m pip install -e packages/$*

install-%: %_contains
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
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
