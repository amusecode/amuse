ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(FORSOCKETS_CFLAGS)
FCFLAGS += -cpp $(DEPFLAGS)
LDFLAGS +=
LDLIBS += $(FORSOCKETS_LIBS)

CODELIB = src/libsphray.a


.PHONY: $(CODELIB)
$(CODELIB):
	$(MAKE) -C src -j 1 amuse_interface FFLAGS="$(FCFLAGS)" FC="$(FC)"


# Building the workers
sphray_worker.f90: interface.py
	amusifier --type=f90 interface.py SPHRayInterface -o $@

sphray_worker.o: sphray_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

sphray_worker: interface.o sphray_worker.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.f90 $(CODELIB)
	$(MPIFC) -o $@ -c $(FCFLAGS) -Isrc $<


# Which packages contain which workers?
amuse-sphray_contains: sphray_worker


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
	$(MAKE) -C src cleanall
	rm -rf *.o *.mod *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
