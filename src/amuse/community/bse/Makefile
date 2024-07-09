ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

# Building the code into a static library
FCFLAGS += $(FORSOCKETS_CFLAGS)
LDLIBS += $(FORSOCKETS_LIBS) 

##### Substitute CODE name #####
BSEOBJ = src/comenv.o src/corerd.o src/deltat.o src/dgcore.o \
         src/evolv2.o src/gntage.o src/hrdiag.o src/instar.o \
         src/kick.o src/mix.o src/mlwind.o src/mrenv.o \
         src/ran3.o src/rl.o src/star.o src/zcnsts.o src/zfuncs.o

##### Modify build command as needed #####
# Building the workers
bse_worker.f90: interface.py
	amusifier --type=f90 interface.py BSEInterface -o $@

bse_worker.o: bse_worker.f90
	$(MPIFC) -c $(FCFLAGS) $<

%.o: %.f
	$(FC) -c -o $@ $(FCFLAGS) $<

bse_worker: interface.o bse_worker.o $(BSEOBJ)
	$(MPIFC) -o bse_worker $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-bse_contains: bse_worker

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
	rm -rf *.o src/*.o *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

