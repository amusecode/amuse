ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
FCFLAGS += $(AMUSE_MPI_CFLAGS) $(FORSOCKETS_CFLAGS)
LDFLAGS +=
LDLIBS += -lm $(AMUSE_MPI_LIBS) $(FORSOCKETS_LIBS)

export FCFLAGS

CODELIB = src/libaarsethzare.a


$(CODELIB): src/*.f
	$(MAKE) -C src all


# Building the workers
cpu_worker.h: interface.py
	amusifier --type=h interface.py AarsethZareInterface -o $@

cpu_worker.f90: interface.py
	amusifier --type=f90 interface.py AarsethZareInterface -o $@

cpu_worker: interface.o cpu_worker.o $(CODELIB)
	$(MPIFC) -o cpu_worker $(LDFLAGS) $^ $(CODELIB) $(LDLIBS)

cpu_worker.o: cpu_worker.f90
	$(MPIFC) -c $(FCFLAGS) $<

%.o: %.f90
	$(FC) -c $(FCFLAGS) $<


# Which packages contain which workers?
amuse-aarsethzare_contains: cpu_worker


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

