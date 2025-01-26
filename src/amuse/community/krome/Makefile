ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
VERSION = 6a762de58a763026e708521ca2b4d68f9a90f669
SHORT_VERSION = $(shell echo $(VERSION) | cut -c 1-12)

krome.tar.gz:
	$(DOWNLOAD) https://bitbucket.org/tgrassi/krome/get/$(VERSION).tar.gz >$@

src/krome: krome.tar.gz
	tar xf $<
	mv tgrassi-krome-$(SHORT_VERSION) src/krome


# Building the code into a static library
DEPFLAGS += $(FORSOCKETS_CFLAGS)
FCFLAGS += $(DEPFLAGS) -Isrc -Isrc/krome/build

LDLIBS += $(FORSOCKETS_LIBS)

CODELIB = src/libchem.a

.PHONY: $(CODELIB)
$(CODELIB): | src/krome
	$(MAKE) -C src libchem.a


# Building the workers
krome_worker.f90: interface.py
	amusifier --type=f90 interface.py KromeInterface -o $@

krome_worker.o: krome_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

krome_worker: interface.o krome_worker.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.f90 $(CODELIB)
	$(MPIFC) -o $@ -c $(FCFLAGS) $<


# Which packages contain which workers?
amuse-krome_contains: krome_worker


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
