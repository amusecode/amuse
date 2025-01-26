ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(FORSOCKETS_CFLAGS)
FCFLAGS += $(DEPFLAGS)
CODE_FFLAGS += -ffixed-form -xf77-cpp-input

LDLIBS += $(FORSOCKETS_LIBS)

CODELIB = src/libnbody6xx.a


$(CODELIB):
	$(MAKE) -C src -j $(CPU_COUNT) amusepar MPIFC="$(MPIFC)" FCFLAGS="$(CODE_FFLAGS) $(FCFLAGS)"


# Building the workers
nbody6xx_worker.f90: interface.py
	amusifier --type=f90 interface.py Nbody6xxInterface -o $@

nbody6xx_worker.o: nbody6xx_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

nbody6xx_worker: interface.o nbody6xx_worker.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.f90
	# interface.f90 includes fixed-form headers from the code, so we need the flags
	$(MPIFC) -o $@ -c $(FCFLAGS) $(CODE_FFLAGS) $<


# Which packages contain which workers?
amuse-nbody6xx_contains: nbody6xx_worker


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
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
