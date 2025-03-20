ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
FCFLAGS += $(FORSOCKETS_CFLAGS)
LDLIBS += $(FORSOCKETS_LIBS) $(LIBS)

MOBSEOBJ = src/comenv.o src/corerd.o src/deltat.o src/dgcore.o \
         src/evolv3.o src/gntage.o src/hrdiag.o src/instar.o \
         src/kick.o src/mix.o src/mlwind.o src/mrenv.o \
         src/ran3.o src/rl.o \
         src/pisn.f src/eddington.f src/fallback.f \
         src/star.o src/zcnsts.o src/zfuncs.o

# Building the workers
mobse_worker.f90: interface.py
	amusifier --type=f90 interface.py MOBSEInterface -o $@

mobse_worker.o: mobse_worker.f90
	$(MPIFC) -c $(FCFLAGS) $<

%.o: %.f
	$(FC) -c -o $@ $(FCFLAGS) $<

mobse_worker: interface.o mobse_worker.o $(MOBSEOBJ)
	$(MPIFC) -o mobse_worker $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-mobse_contains: mobse_worker

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
