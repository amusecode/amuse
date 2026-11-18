ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
VERSION = 8d545530e239b2cceb51c1eafed1e9a8561f8f5f

krome.tar.gz:
	$(DOWNLOAD) https://github.com/tgrassi/krome/archive/$(VERSION).tar.gz >$@

src/krome: krome.tar.gz
	tar xf $<
	mv krome-$(VERSION) src/krome


# Building the code into a static library
DEPFLAGS += $(FORSOCKETS_CFLAGS)
FCFLAGS += $(DEPFLAGS) -Isrc -Isrc/krome/build

LDLIBS += $(FORSOCKETS_LIBS)

CODELIB = src/libchem.a

.PHONY: $(CODELIB)
$(CODELIB): | src/krome
	$(MAKE) -C src libchem.a FC="$(FC)"


# Building the workers
kromesph_worker.f90: interface.py
	amusifier --type=f90 interface.py KromeSphInterface -o $@

kromesph_worker.o: kromesph_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

kromesph_worker: interface.o kromesph_worker.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.f90 $(CODELIB)
	$(MPIFC) -o $@ -c $(FCFLAGS) $<


# Which packages contain which workers?
amuse-kromesph_contains: kromesph_worker


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
	rm -rf *.o *.mod *worker*

.PHONY: distclean
distclean: clean
	$(MAKE) -C src distclean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
