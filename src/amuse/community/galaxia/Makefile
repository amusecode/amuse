ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

# Building the code into a static library
FFLAGS += $(FORSOCKETS_CFLAGS)
LDLIBS += $(FORSOCKETS_LIBS)

LDLIBS += $(LAPACK_LIBS) $(BLAS_LIBS) #-lblas

# Building the workers
interface.o: interface.f90
	$(FC) $(FFLAGS) $(LDLIBS) -c -o $@ $<

galaxia_worker.f90: interface.py
	amusifier --type=f90 interface.py BarAndSpiralsInterface -o $@

galaxia_worker.o: galaxia_worker.f90 interface.o
	$(MPIFC) -c $(FFLAGS) $<

galaxia_worker: galaxia_worker.o interface.o
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

# Which packages contain which workers?
amuse-galaxia_contains: galaxia_worker


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
	rm -rf *.o *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache