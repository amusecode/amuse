ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
##### Remove anything not needed #####
FCFLAGS += $(FORSOCKETS_CFLAGS)
LDLIBS += $(FORSOCKETS_LIBS) $(LIBS)

SSEOBJ = src/deltat.o src/evolv1.o src/hrdiag.o src/kick.o \
		 src/mlwind.o src/mrenv.o src/ran3.o src/star.o \
	  	 src/zcnsts.o src/zfuncs.o

# Building the workers
sse_worker.f90: interface.py
	amusifier --type=f90 interface.py SSEInterface -o $@

sse_worker.o: sse_worker.f90
	$(MPIFC) -c $(FCFLAGS) $<

%.o: %.f
	$(FC) -c -o $@ $(FCFLAGS) $<

sse_worker: interface.o sse_worker.o $(SSEOBJ)
	$(MPIFC) -o sse_worker $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-sse_contains: sse_worker

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
	rm -rf *.o src/*.o *worker* *.mod

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

