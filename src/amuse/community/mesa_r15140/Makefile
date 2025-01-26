ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code and the dependencies

MESA_VERSION := r15140


src:
	mkdir -p src

src/lapack95.tgz: | src
	$(DOWNLOAD) http://user.astro.wisc.edu/~townsend/resource/download/sdk2/src/lapack95.tgz >$@

src/crmath-1.2.zip: | src
	$(DOWNLOAD) https://github.com/rhdtownsend/crmath/archive/v1.2.zip >$@

src/crlibm-1.0beta4.tar.gz: | src
	$(DOWNLOAD) http://user.astro.wisc.edu/~townsend/resource/download/sdk2/src/crlibm-1.0beta4.tar.gz >$@

src/fpx3.tar.gz: | src
	$(DOWNLOAD) http://user.astro.wisc.edu/~townsend/resource/download/sdk2/src/fpx3.tar.gz >$@

src/fpx3_deps: | src
	$(DOWNLOAD) https://raw.githubusercontent.com/rhdtownsend/sdk2/master/profile/common/fpx3/fpx3_deps >$@

src/mesa-$(MESA_VERSION).zip:
	$(DOWNLOAD) https://zenodo.org/record/4311514/files/mesa-$(MESA_VERSION).zip >$@



# Building dependencies into static libraries

## Dependencies are installed into the MESA source directory

MESA_DIR := src/mesa-$(MESA_VERSION)

PATCHES := $(file < patches/series_mesa)
PATCHES := $(patsubst %,%,$(PATCHES))    # replace newlines with spaces


$(MESA_DIR): src/mesa-$(MESA_VERSION).zip
	cd src && unzip -q ../$<
	mkdir -p $(MESA_DIR)/include
	mkdir -p $(MESA_DIR)/lib
	cd $(MESA_DIR) && for p in $(PATCHES) ; do patch -p1 <../../patches/$$p || exit 1 ; done
	# This one should always be skipped
	touch $(MESA_DIR)/adipls/skip_test
	# We can't guarantee bit-for-bit agreement as we don't control the compiler used
	# So everyone has to skip the build tests
	touch $(MESA_DIR)/skip_test


src/LAPACK95: src/lapack95.tgz
	cd src && tar xzf ../$< -m
	cd src && patch -p1 <../patches/lapack95.patch

$(MESA_DIR)/lib/liblapack95.a: | src/LAPACK95 $(MESA_DIR)
	mkdir -p src/LAPACK95/lapack95_modules
	cd src/LAPACK95/SRC && $(MAKE) single_double_complex_dcomplex MPIFC=$(MPIFC)
	cp src/LAPACK95/lapack95_modules/*.mod $(MESA_DIR)/include/
	cp src/LAPACK95/lapack95.a $(MESA_DIR)/lib/liblapack95.a


src/crmath-1.2: src/crmath-1.2.zip
	cd src && unzip -q ../$<

$(MESA_DIR)/lib/libcrmath.a: | src/crmath-1.2
	$(MAKE) -C src/crmath-1.2
	cp src/crmath-1.2/*.mod $(MESA_DIR)/include/
	cp src/crmath-1.2/*.smod $(MESA_DIR)/include/
	cp src/crmath-1.2/libcrmath.a $@


src/crlibm-1.0beta4: src/crlibm-1.0beta4.tar.gz
	cd src && tar xzf ../$<
	cd src && patch -p1 <../patches/crlibm.patch

$(MESA_DIR)/lib/libcrlibm.a: | src/crlibm-1.0beta4
	cd src/crlibm-1.0beta4 && ./configure --enable-static=yes --enable-shared=no
	cd src/crlibm-1.0beta4 && $(MAKE)
	# cd src/crlibm-1.0beta4 && $(MAKE)
	$(AR) cru $@ src/crlibm-1.0beta4/*.o src/crlibm-1.0beta4/scs_lib/*.o


# Building and installing tools

src/bin: | src
	mkdir -p src/bin

src/fpx3: src/fpx3.tar.gz
	cd src && tar xzf ../$< -m

src/bin/fpx3: | src/fpx3 src/bin
	cp src/fpx3/fpx3 src/bin/


src/bin/fpx3_deps: | src/fpx3_deps
	cp $| $@
	cd src/bin && patch <../../patches/fpx3_deps.patch
	chmod u+x $@

## These come with MESA, so no need to download them separately

$(MESA_DIR)/utils/makedepf90-2.8.8: | $(MESA_DIR)
	cd $(MESA_DIR)/utils && tar xzf makedepf90-2.8.8.tar.gz -m

src/bin/makedepf90: | $(MESA_DIR)/utils/makedepf90-2.8.8 src/bin
	cd $(MESA_DIR)/utils/makedepf90-2.8.8 && ./configure
	cd $(MESA_DIR)/utils/makedepf90-2.8.8 && $(MAKE)
	cp $(MESA_DIR)/utils/makedepf90-2.8.8/makedepf90 $@


$(MESA_DIR)/utils/ndiff-2.00: | $(MESA_DIR)
	cd $(MESA_DIR)/utils && tar xzf ndiff-2.00.tar.gz -m
	cp patches/ndiff-configure $(MESA_DIR)/utils/ndiff-2.00/configure

src/bin/ndiff: | $(MESA_DIR)/utils/ndiff-2.00 src/bin
	cd $(MESA_DIR)/utils/ndiff-2.00 && ./configure
	cd $(MESA_DIR)/utils/ndiff-2.00 && $(MAKE)
	cp $(MESA_DIR)/utils/ndiff-2.00/ndiff $@


# Building MESA into a static library

CODELIB := $(MESA_DIR)/lib/libstar.a

MESA_DEPS := $(MESA_DIR)/lib/liblapack95.a
MESA_DEPS += $(MESA_DIR)/lib/libcrmath.a $(MESA_DIR)/lib/libcrlibm.a

MESA_TOOLS := src/bin/fpx3 src/bin/fpx3_deps src/bin/makedepf90 src/bin/ndiff

MESA_ENV := MESA_DIR=$(MESA_DIR) FC=$(FC) MPIFC=$(MPIFC)
MESA_ENV += CC=$(CC) MPICC=$(MPICC) OPENMP_FCFLAGS=$(OPENMP_FCFLAGS)
MESA_ENV += LAPACK_LIBS="$(LAPACK_LIBS)" BLAS_LIBS="$(BLAS_LIBS)"
MESA_ENV += HDF5_FFLAGS="$(HDF5_FFLAGS)" HDF5_LDFLAGS="$(HDF5_LDFLAGS)"

export PATH := $(realpath .)/src/bin:$(PATH)

$(CODELIB): $(MESA_DEPS) $(MESA_TOOLS) | $(MESA_DIR)
	cd $(MESA_DIR) && env $(MESA_ENV) ./install

# Building the workers
DEPFLAGS += $(FORSOCKETS_CFLAGS) $(OPENMP_FFLAGS)
FCFLAGS += $(DEPFLAGS) -I$(MESA_DIR)/include -I$(MESA_DIR)/star/make
FCFLAGS += $(OPENMP_FCFLAGS)

## AMUSE libraries
LDLIBS += $(FORSOCKETS_LIBS)

## MESA libraries
LDFLAGS += $(OPENMP_FCFLAGS) -L$(realpath $(MESA_DIR))/lib $(HDF5_LDFLAGS)
LDLIBS += -lgyre -lionization -latm -lcolors -lstar_data -lnet -leos -lkap -lrates
LDLIBS += -lneu -lchem -linterp_2d -linterp_1d -lnum -lauto_diff -lmtx -lconst
LDLIBS += -lmath -lutils

## MESA dependencies
LDLIBS += -lcrmath -lcrlibm -llapack95
LDLIBS += $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
LDLIBS += -lhdf5_fortran -lhdf5 -lz -ldl


mesa_$(MESA_VERSION)_worker.f90: interface.py
	amusifier --type=f90 interface.py MESAInterface -o $@

mesa_$(MESA_VERSION)_worker.o: mesa_$(MESA_VERSION)_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

mesa_$(MESA_VERSION)_worker: interface.o mesa_interface.o run_star_extras.o mesa_$(MESA_VERSION)_worker.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.f90
	$(MPIFC) -o $@ -c $(FCFLAGS) $<

mesa_interface.o: run_star_extras.o $(CODELIB)

interface.o: mesa_interface.o $(CODELIB)

run_star_extras.o: $(CODELIB)


# Which packages contain which workers?
amuse-mesa_$(MESA_VERSION)_contains: mesa_$(MESA_VERSION)_worker


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
.PHONY: clean_local
clean_local:
	rm -rf *.o *.mod *.smod *worker*

.PHONY: clean
clean: clean_local
	$(MAKE) -C src/LAPACK95/SRC clean
	$(MAKE) -C src/crmath-1.2 clean
	$(MAKE) -C src/crlibm-1.0beta4 clean
	rm -f src/bin/*
	cd src/mesa-$(MESA_VERSION) && ./clean

.PHONY: distclean
distclean: clean_local
	rm -rf src/LAPACK95 src/crmath-1.2 src/crlibm-1.0beta4
	rm -rf src/fpx3 src/bin src/mesa-$(MESA_VERSION)
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
