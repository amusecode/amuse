ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
VERSION = r187

mpiamrvac.tar.gz:
	$(DOWNLOAD) http://amuse.strw.leidenuniv.nl/codes/mpiamrvac-$(VERSION).tgz >$@

PATCHES := $(file < patches/series)
PATCHES := $(patsubst %,%,$(PATCHES))    # replace newlines with spaces

src/mpiamrvac: | mpiamrvac.tar.gz
	tar xf $|
	mkdir -p src
	mv mpiamrvac src/
	for p in $(PATCHES) ; do patch -p1 <patches/$$p ; done



# Building the code into a static library
DEPFLAGS += $(STOPCONDMPI_CFLAGS) $(FORSOCKETS_CFLAGS)

FCFLAGS += $(DEPFLAGS) -std=legacy -ffree-form -fallow-argument-mismatch

LDLIBS += -lm $(STOPCONDMPI_LIBS) $(FORSOCKETS_LIBS)


CODELIB = libmpiamrvac.a

AMRVAC_OPTS_3d := -d=33 -g=14,14,14 -phi=3 -z=2 -p=hd -u=amuse -cp=openmpi
AMRVAC_OPTS_3dacc := -d=33 -g=14,14,14 -phi=3 -z=2 -p=hdacc -u=amuse -cp=openmpi
AMRVAC_OPTS_2d := -d=22 -g=14,14 -phi=3 -z=2 -p=hd -u=amuse -cp=openmpi
AMRVAC_OPTS_2dacc := -d=22 -g=14,14 -phi=3 -z=2 -p=hdacc -u=amuse -cp=openmpi
# should this be -g=14 for a 1D grid? -LV
AMRVAC_OPTS_1d := -d=11 -g=14,14 -phi=3 -z=2 -p=hd -u=amuse -cp=openmpi
AMRVAC_OPTS_1dacc := -d=11 -g=14,14 -phi=3 -z=2 -p=hdacc -u=amuse -cp=openmpi

.PRECIOUS: build%
build%: | src/mpiamrvac
	rm -rf build$*
	cp -r src/mpiamrvac/src build$*
	cd build$* && ./setamrvac $(AMRVAC_OPTS$*)

build%/$(CODELIB): FORCE | build%
	$(MAKE) -C build$* amrvac $(CODELIB) F90FLAGS="$(FCFLAGS)" F90="$(MPIFC)"



# Building the workers
build%/interface.f90: interface.t | build%
	build$*/vacpp.pl $< >$@

build%/interface.o: build%/interface.f90 | build%/$(CODELIB)
	$(MPIFC) -o $@ -c $(FCFLAGS) -Ibuild$* $<


mpiamrvac_worker.f90: interface.py
	amusifier --type=f90 interface.py MpiAmrVacInterface -o $@

.PRECIOUS: build%/mpiamrvac_worker.o
build%/mpiamrvac_worker.o: mpiamrvac_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

mpiamrvac%_worker: build%/interface.o build%/mpiamrvac_worker.o build%/$(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)



# Which packages contain which workers?
WORKERS := mpiamrvac_3d_worker mpiamrvac_3dacc_worker mpiamrvac_2d_worker
WORKERS += mpiamrvac_2dacc_worker mpiamrvac_1d_worker mpiamrvac_1dacc_worker

amuse-mpiamrvac_contains: $(WORKERS)


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
ifneq ($(CI),)
	cd packages/$* && pytest -k 'not noci'
else
	cd packages/$* && pytest
endif



# Cleaning up
.PHONY: clean
clean:
	-$(MAKE) -C src/mpiamrvac/src clean
	rm -f amrvac   # current version of mpiamrvac cleans this too, but it's our dir
	rm -rf build* *.o *.mod *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

FORCE:
