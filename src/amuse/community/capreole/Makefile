ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
FFLAGS_MPI += $(FORSOCKETS_CFLAGS) -fallow-argument-mismatch
FFLAGS += $(FORSOCKETS_CFLAGS) -fallow-argument-mismatch


ifeq  ($(FORTRAN), ifort)
FFLAGS_MPI += -O -ip -u -g -traceback -fpp1 -DIFORT
FFLAGS += -O -ip -u -g -traceback -fpp1 -DIFORT
endif


LDFLAGS += 


LDLIBS_MPI += -lm $(STOPCOND_LIBS) $(FORSOCKETS_LIBS)
LDLIBS += -lm $(STOPCOND_LIBS) $(FORSOCKETS_LIBS)


CODELIB_MPI = src/libcapreole_mpi.a
CODELIB = src/libcapreole.a


$(CODELIB): | src
	$(MAKE) -C src -j 1 amuse_interface F90FLAGS1="$(FFLAGS)" FC="$(FC)"

$(CODELIB_MPI): | src
	$(MAKE) -C src -j 1 amuse_interface_mpi F90FLAGS1="$(FFLAGS_MPI)" FC="$(FC)" MPIFC="$(MPIFC)"


# Building the workers
capreole_worker.f90: interface.py
	amusifier --needs-mpi=FALSE --type=f90 interface.py -i amuse.community.interface.stopping_conditions.StoppingConditionInterface CapreoleInterface -o $@

capreole_worker_mpi.f90: interface.py
	amusifier --type=f90 interface.py -i amuse.community.interface.stopping_conditions.StoppingConditionInterface CapreoleInterface -o $@


capreole_worker.o: capreole_worker.f90
	$(FC) -c $(FFLAGS) $<

capreole_worker_mpi.o: capreole_worker_mpi.f90
	$(MPIFC) -c -o $@ $(FFLAGS_MPI) $<


interface.o: interface.F90 $(CODELIB)
	$(FC) -c $(FFLAGS) -I./src/ -UMPI $<

interface_mpi.o: interface.F90 $(CODELIB_MPI)
	$(MPIFC) -c -o $@ $(FFLAGS_MPI) -I./src/ $<


capreole_worker: interface.o capreole_worker.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

capreole_worker_mpi: interface_mpi.o capreole_worker_mpi.o $(CODELIB_MPI)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS_MPI)


# Which packages contain which workers?
amuse-capreole-mpi_contains: capreole_worker_mpi
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

