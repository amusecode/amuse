ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
SDAR_VERSION = a4ff4b3d076535684313a912ea31985c1431f827
FDPS_VERSION = 55b2bafd316805bad22057cf4ee96217735279bf
PETAR_VERSION = 6ccac364e83ab05e4dad27e6eb6abdf7e5c89bcc

SDAR.tar.gz:
	$(DOWNLOAD) https://github.com/lwang-astro/SDAR/archive/$(SDAR_VERSION).tar.gz >$@

FDPS.tar.gz:
	$(DOWNLOAD) https://github.com/FDPS/FDPS/archive/$(FDPS_VERSION).tar.gz >$@

PeTar.tar.gz:
	$(DOWNLOAD) https://github.com/lwang-astro/PeTar/archive/$(PETAR_VERSION).tar.gz >$@


src/SDAR: SDAR.tar.gz
	mkdir -p src
	cd src && tar xf ../$<
	cd src && mv SDAR-$(SDAR_VERSION) SDAR

src/FDPS: FDPS.tar.gz
	mkdir -p src
	cd src && tar xf ../$<
	cd src && mv FDPS-$(FDPS_VERSION) FDPS

src/PeTar: PeTar.tar.gz
	mkdir -p src
	cd src && tar xf ../$<
	cd src && mv PeTar-$(PETAR_VERSION) PeTar


# Code configuration (it's all header-only, so no CODELIB here)

DEPFLAGS += $(STOPCONDMPI_CFLAGS) $(AMUSE_MPI_CFLAGS)
DEPFLAGS += $(OPENMP_CXXFLAGS)

CXXFLAGS += $(DEPFLAGS)
CXXFLAGS += -Isrc/PeTar/src -Isrc/SDAR/src -Isrc/FDPS/src
CXXFLAGS += -Wall -std=c++17

CXXFLAGS += -D PARTICLE_SIMULATOR_THREAD_PARALLEL
# MPI
CXXFLAGS += -D PARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -D MPICH_IGNORE_CXX_SEEKC

CXXFLAGS += -D SOFT_PERT -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D CLUSTER_VELOCITY
CXXFLAGS += -D USE_QUAD
CXXFLAGS += -D STELLAR_EVOLUTION

CXXFLAGS += -D PROFILE
CXXFLAGS += -D HARD_CHECK_ENERGY
CXXFLAGS += -D TIDAL_TENSOR_3RD

# debug
#CXXFLAGS +=  -D PETAR_DEBUG -D AR_DEBUG -D AR_DEBUG_DUMP -D AR_WARN -D HARD_DEBUG -D HARD_DUMP -D CLUSTER_DEBUG -D ARTIFICIAL_PARTICLE_DEBUG 
CXXFLAGS += -D INTERFACE_DEBUG_PRINT
CXXFLAGS += -D INTERFACE_DEBUG

LDFLAGS += $(OPENMP_CXXFLAGS)

LDLIBS += -lm $(STOPCONDMPI_LIBS) $(AMUSE_MPI_LIBS)
 
src/PeTar/src/get_version.hpp: |src/PeTar
	sed 's/@VERSION@/'`cat src/PeTar/VERSION`'_'`cat src/SDAR/VERSION`'/g' src/PeTar/src/get_version.hpp.in >src/PeTar/src/get_version.hpp


# Building the workers

petar_worker.h: interface.py
	amusifier --type=H interface.py petarInterface -o $@

petar_worker.cc: interface.py
	amusifier --type=c interface.py petarInterface -o $@

petar_worker.o: petar_worker.cc petar_worker.h
	$(MPICXX) -c -o $@ $(CXXFLAGS) $<

petar_worker: petar_worker.o interface.o
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

interface.o: interface.cc src/PeTar/src/get_version.hpp src/PeTar src/SDAR src/FDPS
	$(CXX) -o $@ -c $(CXXFLAGS) $<


# Which packages contain which workers?
amuse-petar_contains: petar_worker


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
	# The tests use up to 5 workers, each of which will use all the cores via OpenMP,
	# slowing things down enormously. By limiting the number of threads, the tests
	# run in a reasonable time at least on my 8C/16T laptop.
	cd packages/$* && OMP_NUM_THREADS=2 pytest


# Cleaning up
.PHONY: clean
clean:
	rm -rf *.o *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
