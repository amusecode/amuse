ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(AMUSE_MPI_CFLAGS)

CFLAGS += $(DEPFLAGS)
FFLAGS += $(DEPFLAGS)

LDFLAGS +=

LDLIBS += -lm $(AMUSE_MPI_LIBS)


CODELIB = src/libsei.a


.PHONY: $(CODELIB)
$(CODELIB): | src
	$(MAKE) -C src -j $(CPU_COUNT) CC='$(MPICC)' all

# Building the workers
sei_worker.h: interface.py
	amusifier --type=h interface.py SeiInterface -o $@

sei_worker.cc: interface.py
	amusifier --type=c interface.py SeiInterface -o $@

interface.o: interface.cc interface.h sei_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) -I./src/ $<

sei_worker.o: sei_worker.cc sei_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

sei_worker: sei_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-sei_contains: sei_worker


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

