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

OMP_FLAGS += $(OPENMP_CFLAGS) $(OPENMP_FFLAGS) $(OPENMP_CXXFLAGS)

OCL_FLAGS += $(CL_CFLAGS)
OCL_FLAGS += -DEVOLVE_OPENCL
INCLUDE += $(CL_FLAGS) -I.

#ifneq (,$(findstring xlc, $(notdir $(CC))))
#CFLAGS += -qlanglvl=extc99
#else
#CFLAGS += -std=gnu99
#endif

LDFLAGS +=

LDLIBS += -lm $(STOPCOND_LIBS) $(STOPCONDMPI_LIBS) $(AMUSE_MPI_LIBS)
LDLIBS += $(FORSOCKETS_LIBS)

OCL_LDLIBS += $(CL_LIBS)


CODELIB = src/libhuayno.a
CODELIB_OPENCL = src/libhuayno_cl.a


$(CODELIB): | src
	$(MAKE) -C src libhuayno.a

$(CODELIB_OPENCL): | src
	$(MAKE) -C src libhuayno_cl.a

interface.o: interface.c
	$(MPICC) -I./src/ -c -o $@ $(CFLAGS) $<

# Building the workers
worker_code.h: interface.py
	amusifier --type=h interface.py HuaynoInterface -o $@

huayno_worker.cc: interface.py
	amusifier --type=c interface.py HuaynoInterface -o $@


huayno_worker.o: huayno_worker.cc worker_code.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

huayno_openmp_worker.o: huayno_worker.cc worker_code.h
	$(MPICXX) -c -o $@ $(CFLAGS) $(OMP_FLAGS) $<

huayno_opencl_worker.o: huayno_worker.cc worker_code.h
	$(MPICXX) -c -o $@ $(CFLAGS) $(OCL_FLAGS) $(INCLUDE) $(OCL_LDLIBS) $<


huayno_worker: huayno_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

huayno_openmp_worker: huayno_openmp_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

huayno_opencl_worker: huayno_opencl_worker.o interface.o $(CODELIB_OPENCL)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)


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

.PHONY: distclean
distclean: clean
	-$(MAKE) -C src distclean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

