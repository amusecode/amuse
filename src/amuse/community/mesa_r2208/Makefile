ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
mesa_r2208.zip:
	$(DOWNLOAD) https://zenodo.org/record/2603017/files/mesa-r2208.zip >$@

PATCHES := $(file < patches/series)
PATCHES := $(patsubst %,%,$(PATCHES))    # replace newlines with spaces

src/mesa: mesa_r2208.zip
	mkdir -p src
	unzip $<
	mv mesa-r2208 src/mesa
	for p in $(PATCHES) ; do patch -p1 <patches/$$p || exit 1 ; done
	cp patches/test_jina.f src/mesa/jina/test/src/test_jina.f


# Building the code into a set of static libraries
DEPFLAGS += $(FORSOCKETS_CFLAGS) $(OPENMP_FFLAGS)

FCFLAGS += $(DEPFLAGS)
LDFLAGS += $(OPENMP_FCFLAGS) -Lsrc/mesa/lib
LDLIBS += $(FORSOCKETS_LIBS)

CODELIB = src/mesa/lib/libstar.a

$(CODELIB): | src/mesa
	cd src/mesa && FC=$(FC) MPIFC=$(MPIFC) ./install

# libstar.a depends on all these other libraries, which are built together with
# it. So we need to link them all.
LDLIBS += -lstar -ldiffusion -lionization -lmlt -latm -lkaro -lcolors
LDLIBS += -leos -lkap -ljina -lweak -lnet -lscreen -lrates -lneu -lchem -lnum
LDLIBS += -linterp_2d -linterp_1d -lutils -lalert -lconst -lmtx -lmesalapack
LDLIBS += -lmesablas

# These test files are modified (see patches/test) and used by the interface (?!)
STAR_TEST = src/mesa/star/test/make

EXTRA_OBJS = $(STAR_TEST)/run_star_support.o $(STAR_TEST)/run_star_extras.o
EXTRA_OBJS += $(STAR_TEST)/calibrate.o $(STAR_TEST)/isochrone.o
EXTRA_OBJS += $(STAR_TEST)/pgstar_stub.o $(STAR_TEST)/create_zams.o
EXTRA_OBJS += $(STAR_TEST)/sample_zams.o $(STAR_TEST)/run_star.o

$(EXTRA_OBJS): $(CODELIB)
# These are built already by the MESA build system when we build libstar.a


# Building the workers
mesa_r2208_worker.f90: interface.py
	amusifier --type=f90 interface.py MESAInterface -o $@

mesa_r2208_worker.o: mesa_r2208_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

mesa_r2208_worker: interface.o mesa_r2208_worker.o $(EXTRA_OBJS)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

interface.o: interface.f90 $(CODELIB)
	$(MPIFC) -o $@ -c $(FCFLAGS) -Isrc/mesa/include -I$(STAR_TEST) -Isrc/mesa/star/make $<

$(STAR_TEST)/%.o:
	@# Keeps make from trying to recompile a .o from the .mod if the .mod is newer.
	@# It thinks the .mod is a Modula-2 program, and the source for the .o...


# Building the data
ZAMS_DATA := data/data/star_data/starting_models/zams_z20m3.data

$(ZAMS_DATA): zams_z20m3.data.gz
	gunzip -c $< >$@

data/data: | src/mesa
	rm -rf data/data
	cp -r src/mesa/data data/


# Which packages contain which workers?
amuse-mesa_r2208_contains: mesa_r2208_worker data/data $(ZAMS_DATA)


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
	cd src/mesa && ./clean
	rm -rf *.o *.mod *worker* data/data $(ZAMS_DATA)

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
