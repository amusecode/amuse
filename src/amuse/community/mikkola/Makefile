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
FFLAGS += $(DEPFLAGS)

LDFLAGS +=

LDLIBS += -lm $(STOPCOND_LIBS) $(STOPCONDMPI_LIBS) $(AMUSE_MPI_LIBS)
LDLIBS += $(FORSOCKETS_LIBS)


# flags from the original MakeFile
#F77FLAGS += 
#ifneq (,$(findstring xlf, $(notdir $(FORTRAN))))
#F77FLAGS += -qfixed -qextname=flush
#endif
#export F77FLAGS


CODELIB = src/libMikkola.a


$(CODELIB): | src
	$(MAKE) -C src all  # make -C src all F77FLAGS="$(F77FLAGS)"


# Building the workers
mikkola_worker.f90: interface.py
	amusifier --type=f90 interface.py MikkolaInterface -o $@

mikkola_worker.o: mikkola_worker.f90
	$(MPIFC) -c $(FFLAGS) $<

interface.o: interface.f90
	$(FC) -c $(FFLAGS) $<

mikkola_worker: mikkola_worker.o interface.o $(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-mikkola_contains: mikkola_worker


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

