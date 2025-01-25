ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
VERSION = ca6907f5f9a95f866b1d7e520cc356ab8cec8dd0

src/phantom.zip:
	mkdir -p src
	$(DOWNLOAD) https://github.com/rieder/phantom/archive/$(VERSION).zip >$@

src/phantom: src/phantom.zip
	unzip $<
	mv phantom-$(VERSION) src/phantom


# Building the code into a static library
FCFLAGS += $(STOPCOND_CFLAGS) $(FORSOCKETS_CFLAGS) $(OPENMP_FCFLAGS)
LDLIBS += -lm $(STOPCOND_LIBS) $(FORSOCKETS_LIBS)
LDFLAGS += $(OPENMP_FCFLAGS)

CODELIB = src/phantom/bin/libphantom-amuse.a

$(CODELIB): | src/phantom
	$(MAKE) -C src/phantom libphantom-amuse SETUP=amuse SYSTEM=gfortran OPENMP=yes


# Building the workers
phantom_worker.f90: interface.py
	amusifier --type=f90 interface.py PhantomInterface -o $@

phantom_worker.o: phantom_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

phantom_worker: interface.o phantom_worker.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.f90
	$(MPIFC) -o $@ -c $(FCFLAGS) $<


# Which packages contain which workers?
amuse-phantom_contains: phantom_worker


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
	$(MAKE) -C src/phantom clean
	rm -rf *.o *.mod *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
