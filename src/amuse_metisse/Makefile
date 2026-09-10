ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
VERSION = 1b79dba5cdb7ac801d4ff24e746a22ad2e566f7c

metisse.tar.gz:
	$(DOWNLOAD) https://github.com/TeamMETISSE/METISSE/archive/$(VERSION).tar.gz >$@

src: metisse.tar.gz
	tar xf $<
	mv metisse-$(VERSION) src

# Building the code into a static library
DEPFLAGS += $(STOPCOND_CFLAGS) $(STOPCONDMPI_CFLAGS) $(AMUSE_MPI_CFLAGS)
DEPFLAGS += $(FORSOCKETS_CFLAGS) $(SIMPLE_HASH_CFLAGS) $(G6LIB_CFLAGS)
DEPFLAGS += $(OPENMP_FFLAGS)
FCFLAGS += $(DEPFLAGS) -ffpe-summary=all -fbacktrace -fcheck=all -g -Wall

LDLIBS += -lm $(STOPCOND_LIBS) $(STOPCONDMPI_LIBS) $(AMUSE_MPI_LIBS)
LDLIBS += $(FORSOCKETS_LIBS)

CODELIB := libmetisse.a

.PHONY: src/$(CODELIB)

src/$(CODELIB): | src/make
	$(MAKE) -C src/make libmetisse


# Building the workers
metisse_worker.f90: interface.py
	amusifier --type=f90 interface.py MetisseInterface -o $@

storage.o: storage.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

metisse_worker.o: metisse_worker.f90 storage.o | src
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

%.o: %.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) -Isrc/make $<

metisse_worker: storage.o src/$(CODELIB) interface.o metisse_worker.o
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-metisse_contains: metisse_worker

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
	-$(MAKE) -C src/make clean
	rm -rf *.o *.mod *worker*

.PHONY: distclean
distclean: clean
	-$(MAKE) -C src/make distclean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

