ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif

# Detecting dependencies
support/config.mk:
	cd support && ./configure


# when tsunami becomes public, replace this block of code
# with a pinned github commit so that the build system
# downloads the archive at install time
tsunami-master.tar.gz:
	@printf "\n"
	@printf "ERROR: '$@' not found in $(CURDIR).\n"
	@printf "\n"
	@printf "  Tsunami is an optional private code and is not distributed with AMUSE.\n"
	@printf "  1. Contact Dr. Alessandro Alberto Trani to obtain the archive.\n"
	@printf "  2. Place it at $(CURDIR)/$@.\n"
	@printf "  3. Re-run the build.\n"
	@printf "\n"
	@exit 1


src/tsunami: tsunami-master.tar.gz
	tar xf $<
	mv tsunami-master src/tsunami

.PHONY: install_tsunami
install_tsunami: | src/tsunami
ifneq ($(VIRTUAL_ENV),)
	pip install "swig>=4.0.2,<4.3" "scikit-build-core>=0.12"
else ifneq ($(CONDA_PREFIX),)
	conda install --yes --freeze-installed --channel conda-forge --override-channels "swig>=4.0.2,<4.3" "scikit-build-core>=0.12"
endif
	cd src/tsunami/ && pip install --no-build-isolation . -Ccmake.define.spin=ON


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
