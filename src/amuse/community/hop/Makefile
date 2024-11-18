ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

# Building the code into a static library
CFLAGS += -DDIFFERENT_MASSES

CODELIB = src/libhop.a

$(CODELIB):
	$(MAKE) -C src -j $(CPU_COUNT) all


# Building the workers
hop_worker.h: interface.py
	amusifier --type=h interface.py HopInterface -o $@

hop_worker.cc: interface.py
	amusifier --type=c interface.py HopInterface -o $@

hop_worker.o: hop_worker.cc hop_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

hop_worker: hop_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

interface.o: interface.cc hop_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<


# Which packages contain which workers?
amuse-hop_contains: hop_worker


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
	rm -rf *.o *worker*

.PHONY: distclean
distclean: clean
	-$(MAKE) -C src clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
