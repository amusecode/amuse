ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
VERSION = c69fa0af018adbbd13d30a7e853c0179b8afbf7f

hermite_grx.tar.gz:
	$(DOWNLOAD) https://github.com/amusecode/Hermite_GRX/archive/$(VERSION).tar.gz >$@

src:
	mkdir -p src

src/Hermite_GRX: | hermite_grx.tar.gz src
	tar xf hermite_grx.tar.gz
	mv Hermite_GRX-$(VERSION) $@


# Building the code into a static library
DEPFLAGS += $(STOPCOND_CFLAGS)
CFLAGS += $(DEPFLAGS) -pthread
LDFLAGS += -pthread

LDLIBS += $(STOPCOND_LIBS)

CODELIB = src/Hermite_GRX/src/libhermite_grx.a


.PHONY: $(CODELIB)
$(CODELIB): | src/Hermite_GRX
	$(MAKE) -C src/Hermite_GRX/src -j $(CPU_COUNT) all CC_FLAGS='$(CFLAGS)' LD_FLAGS='$(LD_FLAGS)'


# Building the workers
hermite_grx_worker.h: interface.py
	amusifier --type=H -i amuse.community.interface.stopping_conditions.StoppingConditionInterface interface.py HermiteGRXInterface -o $@

hermite_grx_worker.cc: interface.py
	amusifier --type=c interface.py HermiteGRXInterface -o $@

hermite_grx_worker.o: hermite_grx_worker.cc hermite_grx_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

hermite_grx_worker: interface.o hermite_grx_worker.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

interface.o: interface.cc | src/Hermite_GRX
	$(MPICXX) -o $@ -c $(CFLAGS) $<


# Which packages contain which workers?
amuse-hermite-grx_contains: hermite_grx_worker


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
	-$(MAKE) -C src/Hermite_GRX/src clean
	rm -rf *.o *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
