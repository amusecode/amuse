ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(STOPCOND_CFLAGS) $(STOPCONDMPI_CFLAGS) $(AMUSE_MPI_CFLAGS)
DEPFLAGS += $(GMP_FLAGS)
DEPFLAGS += $(MPFR_FLAGS)

CFLAGS += $(DEPFLAGS)

# flags from the original Makefile
CXXFLAGS += $(CFLAGS)
CXXFLAGS += -Wall -g -O2  
CXXFLAGS += -std=c++11
CXXFLAGS += -I../mpfrc++


LDFLAGS +=

LDLIBS += -lm $(STOPCOND_LIBS) $(STOPCONDMPI_LIBS) $(AMUSE_MPI_LIBS)

LDLIBS += $(GMP_LIBS)
LDLIBS += $(MPFR_LIBS)


CODELIB = src/libbrutus.a


$(CODELIB): | src
	$(MAKE) -C src -j 1 all CXXFLAGS='$(CXXFLAGS)'


# Building the workers
worker_code.h: interface.py
	amusifier --type=h -i amuse.community.interface.stopping_conditions.StoppingConditionInterface interface.py BrutusInterface -o $@

brutus_worker.cc: interface.py
	amusifier --type=c interface.py BrutusInterface -o $@

interface.o: interface.cc #$(CODELIB)
	$(MPICXX) -c -o $@ $(CFLAGS) -I./src/ -I./mpfrc++/ $<

brutus_worker.o: brutus_worker.cc worker_code.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

brutus_worker: brutus_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-brutus_contains: brutus_worker


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

