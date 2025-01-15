ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

BUILDDIR = build

# Building the code into a static library
DEPFLAGS += $(AMUSE_MPI_CFLAGS)
DEPFLAGS += $(GSL_FLAGS)

CFLAGS += $(DEPFLAGS)

LDFLAGS +=

LDLIBS += -lm $(AMUSE_MPI_LIBS)
LDLIBS += $(GSL_LIBS)

LDLIBS += -L$(BUILDDIR)/src/mmas -l_mmas2 \
	-L$(BUILDDIR)/src/eos -leos \
	-L$(BUILDDIR)/src/usm -lusm \
	-L$(BUILDDIR)/src/std -lstd


$(BUILDDIR)/Makefile:
	@-mkdir $(BUILDDIR)
	cd $(BUILDDIR) && cmake ../src/mmas2 -DGSL_FLAGS="$(GSL_FLAGS)" -DGSL_LIBS="$(GSL_LIBS)"

mmas2: $(BUILDDIR)/Makefile
	$(MAKE) -C $(BUILDDIR) -j $(CPU_COUNT)


# Building the workers
mmams_worker.h: interface.py
	amusifier --type=h interface.py MakeMeAMassiveStarInterface -o $@

mmams_worker.cc: interface.py
	amusifier --type=c interface.py MakeMeAMassiveStarInterface -o $@

interface.o: interface.cc mmams_worker.h mmas2
	$(MPICXX) -c -o $@ $(CFLAGS) -Isrc/mmas2/src $<

mmams_worker.o: mmams_worker.cc mmams_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

mmams_worker: mmams_worker.o interface.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)


# Which packages contain which workers?
amuse-mmams_contains: mmams_worker


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
	rm -rf *.o *worker* $(BUILDDIR)

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

