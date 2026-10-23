ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

# Downloading the code
VERSION = fe9058605a009bebdf6c4b6a107dbcfc57671c96

uclchem.tar.gz:
	$(DOWNLOAD) https://github.com/uclchem/UCLCHEM/archive/$(VERSION).tar.gz >$@

src/UCLCHEM: uclchem.tar.gz
	tar xf $<
	mv UCLCHEM-$(VERSION) src/UCLCHEM


# UCLCHEM is a fortran code which is compiled to a Python API using np.f2py
# For this reason, the UCLCHEM interface only interacts with the Python side
.PHONY: install_uclchem
install_uclchem: | src/UCLCHEM
ifneq ($(CONDA_PREFIX),)
	conda install --yes --channel conda-forge --override-channels pandas numpy pyyaml charset-normalizer matplotlib seaborn scipy tqdm joblib
endif
	pip install src/UCLCHEM


uclchem_worker: interface.py install_uclchem
	PYTHONPATH=. amusifier --type=py --mode=mpi -x interface UclchemInterface UclchemImplementation -o $@

amuse-uclchem_contains: uclchem_worker


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
