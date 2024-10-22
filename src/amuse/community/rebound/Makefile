ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
VERSION = 31d117bdc92182073d0941c331f76e95f515bfc6

$(VERSION).tar.gz:
	$(DOWNLOAD) https://github.com/hannorein/rebound/archive/$(VERSION).tar.gz >$@

src/rebound: | $(VERSION).tar.gz
	tar xf $|
	mv rebound-$(VERSION) src


# Building the code into a static library
CFLAGS += $(OPENMP_CFLAGS) -DCOLLISIONS_DUMMY -DLIBREBOUND
CXXFLAGS += $(STOPCOND_CFLAGS) $(OPENMP_CXXFLAGS)
LDLIBS += $(STOPCOND_LIBS)

CODELIB = src/librebound.a

OBJ_NAMES := rebound tree particle gravity integrator integrator_whfast
OBJ_NAMES += integrator_whfasthelio integrator_ias15 integrator_sei integrator_leapfrog
OBJ_NAMES += integrator_hermes boundary collision tools communication_mpi display
OBJ_NAMES += derivatives glad integrator_janus transformations simulationarchive
OBJ_NAMES += output input

OBJS := $(patsubst %,src/src/%.o,$(OBJ_NAMES))

$(CODELIB): | src/rebound
	$(MAKE) -C src -j $(CPU_COUNT) all CC='$(CC)' CFLAGS='$(CFLAGS)'
	$(AR) rv $@ $(OBJS)
	$(RANLIB) $(CODELIB)


# Building the workers
rebound_worker.h: interface.py
	amusifier --type=h interface.py -i amuse.community.interface.stopping_conditions.StoppingConditionInterface ReboundInterface -o $@

rebound_worker.cc: interface.py
	amusifier --type=c interface.py ReboundInterface -o $@

rebound_worker.o: rebound_worker.cc rebound_worker.h
	$(MPICXX) -c -o $@ $(CXXFLAGS) $(OPENMP_CXXFLAGS) $<

rebound_worker: interface.o rebound_worker.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $(OPENMP_CXXFLAGS) $^ $(LDLIBS)

interface.o: interface.cc rebound_worker.h | src/rebound
	$(MPICXX) -o $@ -c $(CXXFLAGS) -Isrc/src $<


# Which packages contain which workers?
amuse-rebound_contains: rebound_worker


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
	-$(MAKE) -C src clean
	rm -rf *.o *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
