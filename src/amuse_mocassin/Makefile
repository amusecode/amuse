ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
VERSION = 2.02.70

mocassin.$(VERSION).tar.gz:
	$(DOWNLOAD) https://amuse.strw.leidenuniv.nl/codes/mocassin.$(VERSION).tar.gz >$@

PATCHES := $(shell cat patches/series)

src:
	mkdir -p src

src/mocassin: | mocassin.$(VERSION).tar.gz src
	tar xf mocassin.$(VERSION).tar.gz
	mv mocassin.$(VERSION) src/mocassin
	for p in $(PATCHES) ; do patch -p1 <patches/$$p ; done


# Building the code into a static library

# Note that mocassin requires -O2, the tests crash without. Probably some memory
# corruption somewhere...

DEPFLAGS += $(FORSOCKETS_CFLAGS)
FCFLAGS += $(DEPFLAGS) -O2 -std=legacy -fallow-argument-mismatch -fallow-invalid-boz
LDLIBS += $(FORSOCKETS_LIBS)

CODELIB = src/mocassin/libmocassin.a

# Dependencies of interface.f90
MOD_NAMES := common_mod constants_mod dust_mod grid_mod iteration_mod output_mod set_input_mod xsec_mod
MOD_FILES := $(patsubst %,src/mocassin/%.mod,$(MOD_NAMES))


.PHONY: $(CODELIB)
$(CODELIB) $(MOD_FILES) &: | src/mocassin
	$(MAKE) -C src/mocassin -j 1 all FFLAGS='$(FCFLAGS)'


# Building the workers
mocassin_worker.f90: interface.py
	amusifier --type=f90 interface.py MocassinInterface -o $@

mocassin_worker.o: mocassin_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

mocassin_worker: interface.o mocassin_worker.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

interface.o: interface.f90 $(MOD_FILES)
	$(MPIFC) -o $@ -c $(FCFLAGS) -Isrc/mocassin $<


# Which packages contain which workers?
amuse-mocassin_contains: mocassin_worker


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
	$(MAKE) -C src/mocassin clean
	rm -rf *.o *.mod *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
