ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

DEPFLAGS += $(CUDA_FLAGS)

CXXFLAGS += $(DEPFLAGS)

LDFLAGS += $(CUDA_LDFLAGS)


# Building the workers
fastkick_worker.h: interface.py
	amusifier --type=h interface.py FastKickInterface -o $@

fastkick_worker.cc: interface.py
	amusifier --type=c interface.py FastKickInterface -o $@

fastkick_worker.o: fastkick_worker.cc fastkick_worker.h
	$(MPICXX) -c -o $@ $(CXXFLAGS) $<

fastkick_worker: interface.o fastkick_worker.o
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

fastkick_cuda_worker: interface.co cuda_fastkick.co fastkick_worker.o
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS) -lcudart

interface.o interface.co: fastkick_worker.h

%.o: %.cc
	$(MPICXX) -o $@ -c $(CXXFLAGS) $<

%.co: %.cc
	$(MPICXX) -o $@ -c $(CXXFLAGS) -DGPU $<

%.co: %.cu
	$(NVCC) $(NVCCFLAGS) -DGPU -Xcompiler "$(CXXFLAGS)" -c $< -o $@


# Which packages contain which workers?
amuse-fastkick_contains: fastkick_worker

amuse-fastkick-cuda_contains: fastkick_cuda_worker


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
	rm -rf *.o *.co *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
