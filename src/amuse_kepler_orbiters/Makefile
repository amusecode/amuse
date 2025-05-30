ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

# Building the code into a static library
CFLAGS += -Isrc

CODELIB = src/libkepler.a


# Building the code lib
.PHONY: $(CODELIB)
$(CODELIB):
	$(MAKE) -C src -j $(CPU_COUNT) all CXX='$(MPICXX)'

# Building the workers
keplerorbiters_worker.h: interface.py
	amusifier --type=h interface.py KeplerInterface -o $@

keplerorbiters_worker.cc: interface.py
	amusifier --type=c interface.py KeplerInterface -o $@

keplerorbiters_worker.o: keplerorbiters_worker.cc keplerorbiters_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

keplerorbiters_worker: keplerorbiters_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.cc
	$(MPICXX) -c $(CFLAGS) $<


# Which packages contain which workers?
amuse-kepler-orbiters_contains: keplerorbiters_worker


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
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
