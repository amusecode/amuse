ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(STOPCOND_CFLAGS) $(GMP_FLAGS) $(MPFR_FLAGS)

CFLAGS += $(DEPFLAGS)

LDLIBS += $(STOPCOND_LIBS) $(GMP_LIBS) $(MPFR_LIBS)


CODELIB = src/libbrutus.a


$(CODELIB):
	$(MAKE) -C src -j $(CPU_COUNT) all CFLAGS='$(CFLAGS) -I../mpfrc++ -Wall'


# Building the workers
brutus_worker.h: interface.py
	amusifier --type=h interface.py BrutusInterface -o $@

brutus_worker.cc: interface.py
	amusifier --type=c interface.py BrutusInterface -o $@

interface.o: interface.cc # $(CODELIB)
	$(MPICXX) -c -o $@ $(CFLAGS) -Isrc -Impfrc++ $<

brutus_worker.o: brutus_worker.cc brutus_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

brutus_worker: brutus_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-brutus_contains: brutus_worker


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

