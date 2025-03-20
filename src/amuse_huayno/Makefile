ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(STOPCOND_CFLAGS) 
CFLAGS += $(DEPFLAGS)

CFLAGS_CODELIB = $(CFLAGS) -std=gnu99

LDLIBS += -lm $(STOPCOND_LIBS)

CODELIB = libhuayno.a    # relative to build dir
CODELIB_OPENCL = build_opencl/libhuayno_cl.a

build:
	mkdir build

build_%:
	mkdir build_$*

.PHONY: build/$(CODELIB)
build/$(CODELIB): src/Makefile | build
	$(MAKE) -C build -j $(CPU_COUNT) -f ../src/Makefile all VPATH=../src CFLAGS="$(CFLAGS_CODELIB)" CXXFLAGS="$(CXXFLAGS)"

.PHONY: build_openmp/$(CODELIB)
build_openmp/$(CODELIB): src/Makefile | build_openmp
	$(MAKE) -C build_openmp -j $(CPU_COUNT) -f ../src/Makefile all VPATH=../src CFLAGS="$(CFLAGS_CODELIB)" CXXFLAGS="$(CXXFLAGS)" OPENMP_CFLAGS="$(OPENMP_CFLAGS)"

.PHONY: $(CODELIB_OPENCL)
$(CODELIB_OPENCL): src/Makefile | build_opencl
	$(MAKE) -C build_opencl -j $(CPU_COUNT) -f ../src/Makefile libhuayno_cl.a VPATH=../src CFLAGS="$(CFLAGS_CODELIB) $(OPENCL_CFLAGS) -DEVOLVE_OPENCL -I." CXXFLAGS="$(CXXFLAGS)" LIBS="$(OPENCL_LIBS)"


interface.o: interface.c
	$(MPICC) -Isrc -c -o $@ $(CFLAGS) $<


# Building the workers
huayno_worker.h: interface.py
	amusifier --type=h interface.py HuaynoInterface -o $@

huayno_worker.cc: interface.py
	amusifier --type=c interface.py HuaynoInterface -o $@

huayno_worker.o: huayno_worker.cc huayno_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

huayno_openmp_worker.o: huayno_worker.cc huayno_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $(OPENMP_CFLAGS) $<

huayno_opencl_worker.o: huayno_worker.cc huayno_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) -DEVOLVE_OPENCL $(OPENCL_CFLAGS) $<


huayno_worker: huayno_worker.o interface.o build/$(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

huayno_openmp_worker: huayno_openmp_worker.o interface.o build_openmp/$(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $(OPENMP_CFLAGS) $^ $(LDLIBS)

huayno_opencl_worker: huayno_opencl_worker.o interface.o $(CODELIB_OPENCL)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS) $(OPENCL_LIBS)


# Which packages contain which workers?
amuse-huayno_contains: huayno_worker

amuse-huayno-openmp_contains: huayno_openmp_worker

amuse-huayno-opencl_contains: huayno_opencl_worker


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
	rm -rf build*

.PHONY: distclean
distclean: clean
	-$(MAKE) -C src distclean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

