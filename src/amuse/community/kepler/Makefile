ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

# Building the code into a static library
DEPFLAGS += $(AMUSE_MPI_CFLAGS)

CFLAGS += $(DEPFLAGS)
CFLAGS += -Isrc

LDFLAGS +=

LDLIBS += -lm $(AMUSE_MPI_LIBS)

CODELIB = src/libkepler.a


# Buildind the code lib
$(CODELIB): | src
	$(MAKE) -C src all CXX='$(MPICXX)'

# Building the workers
hepler_worker.h: interface.py
	amusifier --type=h interface.py KeplerInterface -o $@

kepler_worker.cc: interface.py
	amusifier --type=c interface.py KeplerInterface -o $@

kepler_worker.o: kepler_worker.cc hepler_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

kepler_worker: kepler_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.cc
	$(MPICXX) -c $(CFLAGS) $<


# Which packages contain which workers?
amuse-kepler_contains: kepler_worker


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

