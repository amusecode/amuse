ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

# Building the code into a static library
FFLAGS += $(FORSOCKETS_CFLAGS)
LDLIBS += $(FORSOCKETS_LIBS)

FFLAGS += -fdefault-real-8 -cpp -fallow-argument-mismatch

CODELIB=src/libfractal.a


# Building the workers
%.o: %.f90
	$(FC) $(FFLAGS) -I./src/ -c -o $@ $<

interface.o: $(CODELIB)

$(CODELIB): src/*.f90
	$(MAKE) -C src amuse_interface

fractalcluster_worker.f90: interface.py
	amusifier --type=f90 interface.py FractalClusterInterface -o $@

fractalcluster_worker.o: fractalcluster_worker.f90
	$(MPIFC) -c $(FFLAGS) $<

fractalcluster_worker: fractalcluster_worker.o interface.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

# Which packages contain which workers?
amuse-fractalcluster_contains: fractalcluster_worker


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

