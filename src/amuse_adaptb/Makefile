ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(STOPCOND_CFLAGS)
DEPFLAGS += $(GMP_FLAGS) $(MPFR_FLAGS)

CFLAGS += $(DEPFLAGS)

LDLIBS += $(STOPCOND_LIBS)

LDLIBS += $(GMP_LIBS)
LDLIBS += $(MPFR_LIBS)


CODELIB = src/libadaptb.a


.PHONY: $(CODELIB)
$(CODELIB):
	$(MAKE) -C src -j $(CPU_COUNT) all


# Building the workers
adaptb_worker.h: interface.py
	amusifier --type=h interface.py AdaptbInterface -o $@

adaptb_worker.cc: interface.py
	amusifier --type=c interface.py AdaptbInterface -o $@

adaptb_worker.o: adaptb_worker.cc adaptb_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

adaptb_worker: interface.o adaptb_worker.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

interface.o: interface.cc adaptb_worker.h
	$(MPICXX) -o $@ -c $(CFLAGS) -Isrc $<


# Which packages contain which workers?
amuse-adaptb_contains: adaptb_worker


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
