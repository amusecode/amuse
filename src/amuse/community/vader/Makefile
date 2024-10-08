ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


VERSION = 0e6357b451e6
VADER_SRC = MJCWilhelm-vader-$(VERSION)

# Downloading the code
vader.tar.gz:
	$(DOWNLOAD) https://bitbucket.org/MJCWilhelm/vader/get/$(VERSION).tar.gz -O $@

# Preparing source directory
src: vader.tar.gz
	mkdir -p src
	tar -xvf $< -C src
	cp -r src/MJCWilhelm-vader-$(VERSION)/vader_csrc/* src
	cp Makefile_vader src/Makefile


# Building the code into a static library
DEPFLAGS += $(AMUSE_MPI_CFLAGS)
DEPFLAGS += $(FORSOCKETS_CFLAGS)
DEPFLAGS += $(GSL_FLAGS)

CFLAGS += $(DEPFLAGS)

LDFLAGS +=

LDLIBS += -lm $(AMUSE_MPI_LIBS)
LDLIBS += $(FORSOCKETS_LIBS)
LDLIBS += $(GSL_LIBS)


CODELIB = src/libvader.a

$(CODELIB): | src
	$(MAKE) -C src -j 1 all MPICC='$(MPICC)' CXXFLAGS='$(CFLAGS)'


# Building the workers
worker_code.h: interface.py
	amusifier --type=h interface.py VaderInterface -o $@

vader_worker.cc: interface.py
	amusifier --type=c interface.py VaderInterface -o $@

interface.o: interface.cc $(CODELIB)
	$(MPICC) -c -o $@ $(CFLAGS) $<

vader_worker.o: vader_worker.cc worker_code.h
	$(MPICC) -c -o $@ $(CFLAGS) $<

vader_worker: vader_worker.o interface.o $(CODELIB)
	$(MPICC) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-vader_contains: vader_worker


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

