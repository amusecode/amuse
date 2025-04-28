ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
CFLAGS += $(STOPCOND_CFLAGS)

LDLIBS += -lm $(STOPCOND_LIBS) -ldl


# Note: there seems to have been a g6/sapporo worker at some point, but it's been
# disabled since at least 2013 and there are no tests for it. Build with -DGPU and see
# the old Makefile if you want to try to make it work again.

CODELIB = build/libbhtree.a


build:
	mkdir build

.PHONY: $(CODELIB)
$(CODELIB): | build
	$(MAKE) -C build -f ../src/Makefile -j $(CPU_COUNT) all VPATH=../src CXXFLAGS="$(CFLAGS)"


# Building the workers
bhtree_worker.h: interface.py
	amusifier --type=h interface.py BHTreeInterface -o $@

bhtree_worker.cc: interface.py
	amusifier --type=c interface.py BHTreeInterface -o $@

bhtree_worker.o: bhtree_worker.cc bhtree_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

bhtree_worker: interface.o bhtree_worker.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

interface.o: interface.cc bhtree_worker.h
	$(MPICXX) -o $@ -c $(CFLAGS) $<


# Which packages contain which workers?
amuse-bhtree_contains: bhtree_worker


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
