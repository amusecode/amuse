ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library

## Flags
DEPFLAGS += $(STOPCOND_CFLAGS) $(FORSOCKETS_CFLAGS)

CFLAGS += $(DEPFLAGS)
FCFLAGS += $(DEPFLAGS) -fdefault-real-8 -ffree-line-length-256 -frecord-marker=4

LDFLAGS +=

LDLIBS += $(STOPCOND_LIBS) $(FORSOCKETS_LIBS) $(FFTW_LIBS) -lm

# The code build uses FFLAGS rather than FCFLAGS, and so does Conda
CODE_MAKE_ARGS = FFLAGS="-I. $(FFLAGS) $(FCFLAGS)" FC="$(FC)" VPATH=../src
CODE_MAKE_ARGS_MP = FFLAGS="-I. $(FFLAGS) $(FCFLAGS) $(OPENMP_FFLAGS)" FC="$(FC)" VPATH="../src"
CODE_MAKE_ARGS_PERIODIC = FFLAGS="-I. $(FFLAGS) $(FCFLAGS) $(FFTW_FLAGS)" FC="$(FC)" VPATH=../src


## Sources and targets
FILIB := libfi.a
MAPLIB := libmap.a


## Building libfi.a plain, with mp and with periodic boundary conditions
build build_mp build_periodic:
	mkdir -p $@

%/globals.h: src/globals_template.h | %
	cp $< $@

build_periodic/globals.h: src/globals_template.h | build_periodic
	cp $< $@
	patch $@ <enable_periodic.patch

MAKE_CODELIB = $(MAKE) -f ../src/Makefile

# We don't know the full dependencies for these targets, so they are phony so that
# we always call src/Makefile for them, which then decides whether to actually do
# anything or not.

.PHONY: build/$(FILIB)
build/$(FILIB): build/globals.h | build 
	$(MAKE_CODELIB) -C build amuse_interface $(CODE_MAKE_ARGS)

.PHONY: build_mp/$(FILIB)
build_mp/$(FILIB): build_mp/globals.h | build_mp
	$(MAKE_CODELIB) -C build_mp amuse_interface $(CODE_MAKE_ARGS_MP)

.PHONY: build_periodic/$(FILIB)
build_periodic/$(FILIB): build_periodic/globals.h | build_periodic
	$(MAKE_CODELIB) -C build_periodic amuse_interface $(CODE_MAKE_ARGS_PERIODIC)


## Building libmap.a with and without mp
.PHONY: build/$(MAPLIB)
build/$(MAPLIB): build/globals.h | build
	$(MAKE_CODELIB) -C build libmap.a $(CODE_MAKE_ARGS)

.PHONY: build_mp/$(MAPLIB)
build_mp/$(MAPLIB): build_mp/globals.h | build_mp
	$(MAKE_CODELIB) -C build_mp libmap.a $(CODE_MAKE_ARGS_MP)


# Building the workers
fi_worker.f90: interface.py
	amusifier --type=f90 --needs-mpi=false $< FiInterface -o $@

fi_map_worker.f90: interface.py
	amusifier --type=f90 --needs-mpi=false $< FiMapInterface -o $@

%.o: %.f90 build/$(FILIB) build/$(MAPLIB)
	$(MPIFC) -c -o $@ $(FCFLAGS) -Ibuild $<

interface_periodic.o: interface.f90 build_periodic/$(FILIB)
	$(MPIFC) -c -o $@ $(FCFLAGS) -Ibuild_periodic $<

fi_worker: interface.o fi_worker.o build/$(FILIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

fi_mp_worker: interface.o fi_worker.o build_mp/$(FILIB)
	$(MPIFC) -o $@ $(LDFLAGS) $(OPENMP_FFLAGS) $^ $(LDLIBS)

fi_map_worker: map_interface.o fi_map_worker.o build/$(MAPLIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

fi_map_mp_worker: map_interface.o fi_map_worker.o build_mp/$(MAPLIB)
	$(MPIFC) -o $@ $(LDFLAGS) $(OPENMP_FFLAGS) $^ $(LDLIBS)

fi_periodic_worker: interface_periodic.o fi_worker.o build_periodic/$(FILIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-fi_contains: fi_worker fi_mp_worker fi_map_worker fi_map_mp_worker fi_periodic_worker

# TODO, if anyone's interested
amuse-fi-gl_contains: fi_gl_worker fi_periodic_gl_worker


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
	$(MAKE) -C src purge
	rm -rf *.o *.mod *worker* build*
	rm -rf tests/output.del tests/__amuse_code_output

.PHONY: distclean
distclean: clean
	$(MAKE) -C src purge
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

