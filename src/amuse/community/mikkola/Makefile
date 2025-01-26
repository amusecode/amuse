ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(STOPCOND_CFLAGS) $(FORSOCKETS_CFLAGS)

FFLAGS += $(DEPFLAGS)

LDLIBS += $(STOPCOND_LIBS) $(FORSOCKETS_LIBS)


CODELIB = src/libMikkola.a


.PHONY: $(CODELIB)
$(CODELIB):
	$(MAKE) -C src -j $(CPU_COUNT) all


# Building the workers
mikkola_worker.f90: interface.py
	amusifier --type=f90 interface.py MikkolaInterface -o $@

mikkola_worker.o: mikkola_worker.f90 interface.o
	$(MPIFC) -c $(FFLAGS) $^

interface.o: interface.f90
	$(FC) -c $(FFLAGS) $<

mikkola_worker: mikkola_worker.o interface.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-mikkola_contains: mikkola_worker


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

