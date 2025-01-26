ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(STOPCOND_CFLAGS)

CFLAGS += $(DEPFLAGS)

LDLIBS += $(STOPCOND_LIBS)

CODELIB = src/libmi6.a
CODELIB_GPU = src/libmi6_gpu.a

.PHONY: $(CODELIB)
$(CODELIB):
	$(MAKE) -C src -j $(CPU_COUNT) emulib

.PHONY: $(CODELIB_GPU)
$(CODELIB_GPU):
	$(MAKE) -C src -j $(CPU_COUNT) gpulib


# Building the workers
mi6_worker.h: interface.py
	amusifier --type=h interface.py MI6Interface -o $@

mi6_worker.cc: interface.py
	amusifier --type=c interface.py MI6Interface -o $@

mi6_worker.o: mi6_worker.cc mi6_worker.h
	$(MPICXX) -c -o $@ -Isrc $(CFLAGS) $<

mi6_worker: interface.o mi6_worker.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

# TODO: Re-enable and make work when we have Sapporo2 packaged and available
mi6_sapporo2_worker: mi6_worker.o interface.o $(CODELIB_GPU)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS) $(SAPPORO2_LIBS)

interface.o: interface.cc mi6_worker.h
	$(MPICXX) -o $@ -c $(CFLAGS) -Isrc $<


# Which packages contain which workers?
amuse-mi6_contains: mi6_worker

amuse-mi6-sapporo2_contains: mi6_sapporo2_worker


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
	rm -rf *.o *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
