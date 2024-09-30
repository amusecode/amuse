ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

# Building the code into a static library
DEPFLAGS += $(OPENMP_CXXFLAGS) $(CUDA_FLAGS)
CXXFLAGS += -O3 -Wall -DCHECK_ERRORS $(DEPFLAGS)

# Custom options, should really be run-time
# CXXFLAGS += -DCHECK_TIMES -DPLUMMER -DUNROLL -DGPUCORR -DGALAXY

NVCCFLAGS := -O3 -Xcompiler "$(CXXFLAGS) $(MPI_CXXFLAGS)" -Ilib

LDFLAGS += $(OPENMP_CXXFLAGS) $(CUDA_LDFLAGS)

LDLIBS += -lcudart

CODELIB = src/libhigpus.a

$(CODELIB):
	$(MAKE) -C src -j $(CPU_COUNT) all NVCC="$(NVCC)" NVCCFLAGS='$(NVCCFLAGS)' CXX="$(MPICXX)"


# Building the workers
higpus_cuda_worker.h: interface.py
	amusifier --type=H interface.py HiGPUsInterface -o $@

higpus_cuda_worker.cc: interface.py
	amusifier --type=c interface.py HiGPUsInterface -o $@

higpus_cuda_worker.o: higpus_cuda_worker.cc higpus_cuda_worker.h
	$(MPICXX) -c -o $@ $(CXXFLAGS) $<

higpus_cuda_worker: higpus_cuda_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.cc
	$(MPICXX) -o $@ -c $(CXXFLAGS) $<


# Which packages contain which workers?
amuse-higpus-cuda_contains: higpus_cuda_worker


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
