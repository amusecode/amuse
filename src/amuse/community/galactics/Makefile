ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a set of binaries
CFLAGS += $(DEPFLAGS)
FCFLAGS += -ffixed-line-length-0 -fno-backslash -fallow-argument-mismatch

EXENAMES := dbh genhalo genbulge getfreqs gendisk diskdf toascii
GASEXENAMES := $(EXENAMES) gengas

EXES := $(patsubst %,data/bin/%,$(EXENAMES))
GASEXES := $(patsubst %,data/gbin/%,$(GASEXENAMES))

data/bin:
	mkdir -p data/bin

data/gbin:
	mkdir -p data/gbin

$(EXES) &: | data/bin
	$(MAKE) -C src/src all CC="$(CC)" CFLAGS="$(CFLAGS)" F77="$(FC)" FFLAGS="$(FCFLAGS)" PGSTUB=1
	$(MAKE)  -C src/src install
	$(MAKE) -C src/potsrc all CC="$(CC)" CFLAGS="$(CFLAGS)" F77="$(FC)" FFLAGS="$(FCFLAGS)" PGSTUB=1
	$(MAKE) -C src/potsrc install
	cp src/bin/* data/bin/

$(GASEXES) &: | data/gbin
	$(MAKE) -C gas_src/src all CC="$(CC)" CFLAGS="$(CFLAGS)" F77="$(FC)" FFLAGS="$(FCFLAGS)" FFLAGS2="$(FCFLAGS)" PGSTUB=1
	$(MAKE) -C gas_src/src install
	cp gas_src/bin/* data/gbin/


# Which packages contain which (pseudo-)workers?
amuse-galactics_contains: $(EXES) $(GASEXES)


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
	$(MAKE) -C src/src clean
	$(MAKE) -C src/potsrc clean
	$(MAKE) -C gas_src/src clean
	rm -rf packages/amuse-galactics/interface packages/amuse-galactics/__amuse_code_output
	rm -rf *worker* data

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
