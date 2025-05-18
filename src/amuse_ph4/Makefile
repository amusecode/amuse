ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
CFLAGS += $(STOPCOND_CFLAGS)
LDLIBS += -lm $(STOPCOND_LIBS)

CODELIB = src/libph4.a
CODELIB_SAPPORO = src/libph4gpu.a


.PHONY: $(CODELIB)
$(CODELIB):
	$(MAKE) -C src -j $(CPU_COUNT) all CXX='$(MPICXX)' CFLAGS='$(CFLAGS)'

.PHONY: $(CODELIB_SAPPORO)
$(CODELIB_SAPPORO):
	$(MAKE) -C src -j $(CPU_COUNT) gpulib CXX='$(MPICXX)' CFLAGS='$(CFLAGS)' SAPPORO_LIBDIRS='$(SAPPORO_LIGHT_CFLAGS)' SAPPORO_LIBS='$(SAPPORO_LIGHT_LIBS)'


# Building the workers

ph4_worker.h: interface.py
	amusifier --type=h interface.py ph4Interface -o $@

ph4_worker.cc: interface.py
	amusifier --type=c interface.py ph4Interface -o $@

ph4_worker.o: ph4_worker.cc ph4_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

ph4_worker: ph4_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

ph4_sapporo_worker: ph4_worker.o interface.o $(CODELIB_SAPPORO)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS) $(SAPPORO_LIGHT_LIBS)

%.o: %.cc
	$(MPICXX) -o $@ -c $(CFLAGS) $<


# Which packages contain which workers?
amuse-ph4_contains: ph4_worker

amuse-ph4-sapporo_contains: ph4_sapporo_worker


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
	$(MAKE) -C src clean
	rm -rf *.o *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

