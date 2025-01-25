ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
FFLAGS += -fallow-argument-mismatch $(FORSOCKETS_CFLAGS)
FCFLAGS = $(FFLAGS)

export CFLAGS
export FFLAGS
export LDFLAGS

FCFLAGS_NOMPI := $(FCFLAGS) $(G6LIB_CFLAGS) $(STOPCONDMPI_CFLAGS) -DNOMPI
FCFLAGS_MPI := $(FCFLAGS) $(G6LIB_CFLAGS) $(STOPCONDMPI_CFLAGS)
FCFLAGS_SAPPORO := $(FCFLAGS) $(SAPPORO_LIGHT_CFLAGS) $(STOPCOND_CFLAGS) -DNOMPI

LDLIBS += $(FORSOCKETS_LIBS)

LDLIBS_NOMPI := $(LDLIBS) $(STOPCONDMPI_LIBS) $(G6LIB_LIBS)
LDLIBS_MPI := $(LDLIBS) $(STOPCONDMPI_LIBS) $(G6LIB_LIBS)
LDLIBS_SAPPORO += $(LDLIBS) $(STOPCOND_LIBS) $(SAPPORO_LIGHT_LIBS)

CODELIB_NOMPI = build_nompi/libphigrape.a
CODELIB_MPI = build_mpi/libphigrape.a
CODELIB_SAPPORO = build_sapporo/libphigrape.a


build_nompi:
	mkdir -p $@

$(CODELIB_NOMPI): | build_nompi
	$(MAKE) -C build_nompi -f ../src/Makefile -j 1 all FC="$(FC)" CC="$(CC)" CFLAGS="$(CFLAGS) -DNOMPI" FFLAGS="$(FFLAGS) -DNOMPI" OPTS="-DNOMPI" STOPCOND_CFLAGS="$(STOPCOND_CFLAGS)" STOPCOND_LIBS="$(STOPCOND_LIBS)" G6LIBS="$(G6LIB_LIBS)" VPATH=../src

build_mpi:
	mkdir -p $@

$(CODELIB_MPI): | build_mpi
	$(MAKE) -C build_mpi -f ../src/Makefile -j 1 all FC="$(MPIFC)" CC="$(MPICC)" CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)" STOPCOND_CFLAGS="$(STOPCONDMPI_CFLAGS)" STOPCOND_LIBS="$(STOPCONDMPI_LIBS)" G6LIBS="$(G6LIB_LIBS)" VPATH=../src

build_sapporo:
	mkdir -p $@

$(CODELIB_SAPPORO): | build_sapporo
	$(MAKE) -C build_sapporo -f ../src/Makefile -j 1 all FC="$(FC)" CC="$(CC)" CFLAGS="$(CFLAGS) -DNOMPI" FFLAGS="$(FFLAGS) -DNOMPI"  STOPCOND_CFLAGS="$(STOPCOND_CFLAGS)" STOPCOND_LIBS="$(STOPCOND_LIBS)" VPATH=../src


# Building the workers
phigrape_worker.f90: interface.py
	amusifier --type=f90 interface.py PhiGRAPEInterface -o $@


phigrape_worker.o: phigrape_worker.f90
	$(MPIFC) -c -o $@ -Isrc $(FCFLAGS) $<

phigrape_worker: interface.o phigrape_worker.o $(CODELIB_NOMPI)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS_NOMPI)

interface.o: interface.F
	$(MPIFC) -o $@ -c -Isrc $(FCFLAGS_NOMPI) $<


phigrape_mpi_worker.o: phigrape_worker.f90
	$(MPIFC) -c -o $@ -Isrc $(FCFLAGS_MPI) $<

phigrape_mpi_worker: phigrape_mpi_worker.o interface_mpi.o $(CODELIB_MPI)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS_MPI)

interface_mpi.o: interface.F
	$(MPIFC) -o $@ -c -Isrc $(FCFLAGS_MPI) $<


phigrape_sapporo_worker.o: phigrape_worker.f90
	$(MPIFC) -c -o $@ -Isrc $(FCFLAGS_SAPPORO) $<

phigrape_sapporo_worker: interface_sapporo.o phigrape_sapporo_worker.o $(CODELIB_SAPPORO)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS_SAPPORO)

interface_sapporo.o: interface.F
	$(MPIFC) -o $@ -c -Isrc $(FCFLAGS_SAPPORO) $<


# Which packages contain which workers?
amuse-phigrape_contains: phigrape_worker phigrape_mpi_worker

amuse-phigrape-sapporo_contains: phigrape_sapporo_worker


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
	rm -rf build_nompi build_mpi build_sapporo

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
