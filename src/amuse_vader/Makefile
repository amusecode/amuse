ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


VERSION = 0e6357b451e6

# Downloading the code
vader.tar.gz:
	$(DOWNLOAD) https://bitbucket.org/MJCWilhelm/vader/get/$(VERSION).tar.gz >$@

# Preparing source directory
src: vader.tar.gz
	mkdir -p src
	tar -xvf $< -C src
	cp -r src/MJCWilhelm-vader-$(VERSION)/vader_csrc/* src
	rm -f src/Makefile
	ln -s ../Makefile_vader src/Makefile


# Building the code into a static library
PROBS := none cmzdisk cmzdisk2 gidisk pedisk ring ringrad selfsim
USERFUNCS := $(patsubst %,src/prob/userFunc_%.o,$(PROBS))

DEPFLAGS += $(GSL_FLAGS)

CFLAGS += $(DEPFLAGS)

LDFLAGS +=

LDLIBS += $(GSL_LIBS)


CODELIB = src/libvader.a

.PHONY: $(CODELIB)
$(CODELIB) $(USERFUNCS) &: | src
	$(MAKE) -C src -j $(CPU_COUNT) all MPICC='$(MPICC)' CXXFLAGS='$(CFLAGS)' PROBS='$(PROBS)'


# Building the workers

WORKERS := $(patsubst %,vader_%_worker,$(PROBS))

vader_worker.h: interface.py
	amusifier --type=H interface.py VaderInterface -o $@

vader_worker.cc: interface.py
	amusifier --type=c interface.py VaderInterface -o $@

interface.o: interface.cc $(CODELIB)
	$(MPICXX) -c -o $@ $(CFLAGS) $<

vader_worker.o: vader_worker.cc vader_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

vader_%_worker: vader_worker.o interface.o $(CODELIB) src/prob/userFunc_%.o
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-vader_contains: $(WORKERS)


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

