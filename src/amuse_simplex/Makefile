ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(GSL_FLAGS) $(HDF5_CFLAGS) $(HEALPIX_CXX_CFLAGS) $(QHULL_FLAGS)

CFLAGS += $(DEPFLAGS)

LDFLAGS += $(HDF5_LDFLAGS)

LDLIBS += $(GSL_LIBS) $(HDF5_LIBS) $(HEALPIX_CXX_LIBS) $(QHULL_LIBS)

CODELIB = libsimplex.a

CODEOBJS = amuse_interface.o SimpleX.o Structs.o rates.o Common.o h5w_serial.o
CODEOBJS += keyValue.o tree_structures.o hilbert.o
CODEOBJS := $(patsubst %,src/obj/%,$(CODEOBJS))

.PHONY: $(CODELIB)
$(CODELIB):
	$(MAKE) -C src amuse MPICXX=$(MPICXX) HDF5_FLAGS="$(HDF5_CFLAGS) $(HDF5_CPPFLAGS)" HDF5_LIBS="$(HDF5_LDFLAGS) $(HDF5_LIBS)"
	$(AR) ruv $@ $(CODEOBJS)
	$(RANLIB) $@


# Building the workers
simplex_worker.h: interface.py
	amusifier --type=h interface.py SimpleXInterface -o $@

simplex_worker.cc: interface.py
	amusifier --type=c interface.py SimpleXInterface -o $@

simplex_worker.o: simplex_worker.cc simplex_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

simplex_worker: simplex_worker.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-simplex_contains: simplex_worker


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
