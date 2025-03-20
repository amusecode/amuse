ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the workers

# mameclot is a stand-alone program that is run as a subprocess from interface.py,
# communicating via pipes. So no amusifier here and no static library either.

CFLAGS += $(CUDA_FLAGS)
LDFLAGS += $(CUDA_LDFLAGS)

mameclot_worker:
	$(MAKE) -C src -j $(CPU_COUNT) CC="$(CC)" CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)"
	cp src/mameclot mameclot_worker


# Which packages contain which workers?
amuse-mameclot_contains: mameclot_worker

# There is some CUDA code here for a mameclot_cuda_worker, but it's not used by
# interface.py. For the future. - LV

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
	$(MAKE) -C src clean
	rm -rf *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
