ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

# Building the code into a static library
DEPFLAGS += $(STOPCOND_CFLAGS)

CFLAGS += $(DEPFLAGS)
# CFLAGS += -Wall -g -O2

LDFLAGS +=

LDLIBS += -lm $(STOPCOND_LIBS)

CODELIB = src/libsakura.a


.PHONY: $(CODELIB)
$(CODELIB): | src
	$(MAKE) -C src -j 1 CXX='$(MPICXX)' CXXFLAGS='$(CFLAGS)' all

interface.o: interface.cc sakura_worker.h
	$(MPICXX) -I./src/ -c -o $@ $(CFLAGS) $<

# Building the workers
sakura_worker.h: interface.py
	amusifier --type=h -i amuse.community.interface.stopping_conditions.StoppingConditionInterface interface.py SakuraInterface -o $@

sakura_worker.cc: interface.py
	amusifier --type=c interface.py SakuraInterface -o $@

sakura_worker.o: sakura_worker.cc sakura_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

sakura_worker: sakura_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-sakura_contains: sakura_worker


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

