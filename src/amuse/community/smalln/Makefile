ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(STOPCOND_CFLAGS) $(STOPCONDMPI_CFLAGS) $(AMUSE_MPI_CFLAGS)
DEPFLAGS += $(FORSOCKETS_CFLAGS)

CFLAGS += $(DEPFLAGS)
CFLAGS += -Wall -g -Isrc

LDFLAGS +=

LDLIBS += -lm $(STOPCOND_LIBS) $(STOPCONDMPI_LIBS) $(AMUSE_MPI_LIBS)
LDLIBS += $(FORSOCKETS_LIBS)


CODELIB = src/libsmallN.a


$(CODELIB): | src
	$(MAKE) -C src all CXX='$(MPICXX)'

# Building the workers
smalln_worker.h: interface.py
	amusifier --type=h interface.py -i amuse.community.interface.stopping_conditions.StoppingConditionInterface SmallNInterface -o $@

smalln_worker.cc: interface.py
	amusifier --type=c interface.py SmallNInterface -o $@

interface.o: interface.cc
	$(MPICXX) -c -o $@ $(CFLAGS) $<

smalln_worker.o: smalln_worker.cc smalln_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

smalln_worker: smalln_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-smalln_contains: smalln_worker


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

