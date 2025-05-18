ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
FCFLAGS += $(FORSOCKETS_CFLAGS) -fallow-argument-mismatch

LDLIBS += $(STOPCOND_LIBS) $(FORSOCKETS_LIBS)


CODELIB = src/libcapreole_mpi.a

.PHONY: $(CODELIB)
$(CODELIB):
	$(MAKE) -C src -j 1 amuse_interface_mpi F90FLAGS1="$(FCFLAGS)" FC="$(FC)"


# Building the workers
capreole_worker.f90: interface.py
	amusifier --type=f90 interface.py -i amuse.community.interface.stopping_conditions.StoppingConditionInterface CapreoleInterface -o $@

capreole_worker.o: capreole_worker.f90
	$(MPIFC) -c -o$@ $(FCFLAGS) $<

interface.o: interface.F90 $(CODELIB)
	$(FC) -c -o $@ $(FCFLAGS) -I./src/ -DMPI $<

capreole_worker: interface.o capreole_worker.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-capreole_contains: capreole_worker


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
ifneq ($(CI),)
	cd packages/$* && pytest -k noci
else
	cd packages/$* && pytest
endif



# Cleaning up
.PHONY: clean
clean:
	-$(MAKE) -C src clean
	rm -rf *.o *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

