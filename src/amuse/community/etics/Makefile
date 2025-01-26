ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(CUDA_FLAGS)
CXXFLAGS += $(DEPFLAGS) -Wall
LDFLAGS += $(CUDA_LDFLAGS)
LDLIBS += -lcudart

CODELIB = src/libetics.a

src/config.mk:
	cd src && ../src/configure

.PHONY: $(CODELIB)
$(CODELIB): src/config.mk
	$(MAKE) -C src -j $(CPU_COUNT) libetics.a


# Building the workers
etics_worker.h: interface.py
	amusifier --type=H interface.py EticsInterface -o $@

etics_worker.cc: interface.py
	amusifier --type=c interface.py EticsInterface -o $@

etics_worker.o: etics_worker.cc etics_worker.h
	$(MPICXX) -c -o $@ $(CXXFLAGS) $<

etics_cuda_worker: interface.o etics_worker.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.cu
	$(NVCC) -o $@ -c -Xcompiler="$(CXXFLAGS)" $<

# Which packages contain which workers?
amuse-etics-cuda_contains: etics_cuda_worker


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
