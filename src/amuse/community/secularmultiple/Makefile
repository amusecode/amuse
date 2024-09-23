ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(AMUSE_MPI_CFLAGS)
DEPFLAGS += $(FORSOCKETS_CFLAGS)

CFLAGS += $(DEPFLAGS)

CFLAGS += -g -DTOOLBOX -O3

LDFLAGS +=

LDLIBS += -lm $(AMUSE_MPI_LIBS)
LDLIBS += $(FORSOCKETS_LIBS)

SCR = src/

OBJS = src/types.o src/evolve.o src/structure.o src/ODE_system.o src/root_finding.o src/newtonian.o src/postnewtonian.o src/tides.o src/external.o src/cvode/cvode.o src/cvode/cvode_dense.o src/cvode/cvode_direct.o src/cvode/cvode_io.o src/cvode/nvector_serial.o src/cvode/sundials_dense.o src/cvode/sundials_direct.o src/cvode/sundials_math.o src/cvode/sundials_nvector.o 

# Building the workers
secularmultiple_worker.h: interface.py
	amusifier --type=h interface.py SecularMultipleInterface -o $@

secularmultiple_worker.cc: interface.py
	amusifier --type=c interface.py SecularMultipleInterface -o $@

interface.o: interface.cpp
	$(MPICXX) -I.$(SCR) -c -o $@ $(CFLAGS) $<

secularmultiple_worker.o: secularmultiple_worker.cc
	$(MPICXX) -I.$(SCR) -c -o $@ $(CFLAGS) $<

secularmultiple_worker: secularmultiple_worker.o interface.o $(OBJS)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.cpp
	$(MPICXX) -c -o $@ $(CFLAGS) $<


# Which packages contain which workers?
amuse-secularmultiple_contains: secularmultiple_worker


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
	rm -rf src/*.o

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

