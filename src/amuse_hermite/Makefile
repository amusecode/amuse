ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

# Building the code into a static library
DEPFLAGS += $(STOPCONDMPI_CFLAGS) $(AMUSE_MPI_CFLAGS)

CFLAGS += $(DEPFLAGS)

LDFLAGS +=

LDLIBS += -lm $(STOPCONDMPI_LIBS) $(AMUSE_MPI_LIBS)


# Building the workers
hermite_worker.h: interface.py
	amusifier --type=h -i amuse.community.interface.stopping_conditions.StoppingConditionInterface interface.py HermiteInterface -o $@

hermite_worker.cc: interface.py
	amusifier --type=c interface.py HermiteInterface -o $@

hermite_worker.o: hermite_worker.cc hermite_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

hermite_worker: hermite_worker.o interface.o
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.cc
	$(MPICXX) -c $(CFLAGS) -o $@ $<

# Which packages contain which workers?
amuse-hermite_contains: hermite_worker


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
	cd packages/$* && pytest -k 'not noci'
else
	cd packages/$* && pytest
endif



# Cleaning up
.PHONY: clean
clean:
	rm -rf *.o *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

