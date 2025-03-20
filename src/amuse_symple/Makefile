ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(STOPCONDMPI_CFLAGS) $(AMUSE_MPI_CFLAGS)

CFLAGS += $(DEPFLAGS)
CFLAGS += -O3 -Wall -DTOOLBOX

LDFLAGS +=

LDLIBS += -lm $(STOPCONDMPI_LIBS) $(AMUSE_MPI_LIBS)


# Building the workers
symple_worker.h: interface.py
	amusifier --type=h -i amuse.community.interface.stopping_conditions.StoppingConditionInterface interface.py sympleInterface -o $@

symple_worker.cc: interface.py
	amusifier --type=c interface.py sympleInterface -o $@

interface.o: interface.cc symple_worker.h vec3.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

symple_worker.o: symple_worker.cc symple_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

symple_worker: symple_worker.o interface.o
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-symple_contains: symple_worker


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
	rm -rf *.o *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

