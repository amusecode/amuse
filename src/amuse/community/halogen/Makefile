ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

# Building the code into a static library
LDLIBS += -lm

CODELIB = src/libhalogen.a

$(CODELIB):
	$(MAKE) -C src -j $(CPU_COUNT) all CC='$(CC)'

interface.o: interface.cc halogen_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

# Building the workers
halogen_worker.h: interface.py
	amusifier --type=h interface.py HalogenInterface -o $@

halogen_worker.cc: interface.py
	amusifier --type=c interface.py HalogenInterface -o $@

halogen_worker.o: halogen_worker.cc halogen_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

halogen_worker: halogen_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-halogen_contains: halogen_worker


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
	-$(MAKE) -C src clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

