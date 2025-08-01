ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(STOPCOND_CFLAGS) $(FORSOCKETS_CFLAGS)
FCFLAGS += $(DEPFLAGS)

LDLIBS += $(STOPCOND_LIBS) $(FORSOCKETS_LIBS)

CODELIB = src/libmercury.a


.PHONY: $(CODELIB)
$(CODELIB):
	$(MAKE) -C src -j $(CPU_COUNT) amuse_interface FFLAGS="$(FCFLAGS)" FC="$(FC)" F77FLAGS="$(F77FLAGS)"


# Building the workers
mercury_worker.f90: interface.py
	amusifier --type=f90 interface.py MercuryInterface -o $@

mercury_worker.o: mercury_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

mercury_worker: interface.o mercury_worker.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.f90 | $(CODELIB)
	$(MPIFC) -o $@ -c $(FCFLAGS) -Isrc $<


# Which packages contain which workers?
amuse-mercury_contains: mercury_worker


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
	rm -rf *.o *.mod *worker* $(CODELIB)

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
