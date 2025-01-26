ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

# Building the code into static libraries
DEPFLAGS += $(STOPCOND_CFLAGS) $(AMUSE_MPI_CFLAGS) $(FFTW_FLAGS)

CFLAGS += $(DEPFLAGS) -DAMUSE

LDLIBS += -lm $(STOPCOND_LIBS) $(AMUSE_MPI_LIBS) $(FFTW_LIBS)

CODELIB = bin/libathena.a    # relative to build dir

CONFIGOPTS_normal = --enable-mpi --with-gas=hydro --enable-smr
CONFIGOPTS_selfgrav = --enable-mpi --with-gas=hydro --with-gravity=fft --enable-fft --disable-cooling
CONFIGOPTS_mhd = --enable-mpi --with-gas=mhd --with-eos=adiabatic --enable-smr --with-order=3 --with-flux=hlld
CONFIGOPTS_scalar = --enable-mpi --with-gas=hydro --enable-smr --with-nscalars=1

build_%:
	cp -r src/athena build_$*
	cd build_$* && ln -fs ../../amuse_problem_4.1.c src/problem.c

build_%/Makefile: | build_%
	cd build_$* && ./configure $(CONFIGOPTS_$*)

.PHONY: build_normal/$(CODELIB) build_selfgrav/$(CODELIB) build_mhd/$(CODELIB) build_scalar/$(CODELIB)
build%/$(CODELIB): | build%/Makefile
	$(MAKE) -C build$* -j $(CPU_COUNT) all CC=$(MPICC) LDR=$(MPICC) OPT="$(CFLAGS)" FFTW_FLAGS="$(FFTW_FLAGS)" FFTW_LIBS="$(FFTW_LIBS)"


# Building the workers
athena_worker.h: interface.py
	amusifier --type=h interface.py AthenaInterface -o $@

athena_worker.cc: interface.py
	amusifier --type=c interface.py AthenaInterface -o $@

athena_worker.o: athena_worker.cc athena_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

athena_%_worker: build_%/interface_4.1.o athena_worker.o build_%/$(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

build_%/interface_4.1.o: interface_4.1.c athena_worker.h | build_%/Makefile
	$(MPICC) -o $@ -c $(CFLAGS) -Ibuild_$*/src $<


# Which packages contain which workers?
amuse-athena_contains: athena_normal_worker athena_selfgrav_worker athena_mhd_worker athena_scalar_worker


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
	-$(MAKE) -C src/athena/src clean
	rm -rf build*
	rm -rf *.o *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
