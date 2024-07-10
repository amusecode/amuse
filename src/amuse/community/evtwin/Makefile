ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure

# Building the code into a static library
FFLAGS += $(FORSOCKETS_CFLAGS) -cpp
LDFLAGS +=
LDLIBS += $(FORSOCKETS_LIBS)

SRCDIR = src/trunk
EVTWINLIB=$(SRCDIR)/lib/libtwin.a $(SRCDIR)/lib/libtwincore.a
BUILDDIR = build
CMAKE_DEPS = $(SRCDIR)/CMakeLists.txt $(SRCDIR)/cmake/CompilerFlags.cmake $(SRCDIR)/cmake/FileList.cmake
CMAKE_GENERATOR ?= Unix Makefiles


$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/Makefile: $(CMAKE_DEPS) | $(BUILDDIR)
	cd $(BUILDDIR) && FC="$(FC) $(FCFLAGS)" cmake ../src/trunk -DCMAKE_Fortran_FLAGS:STRING="$(FCFLAGS)"  -G "$(CMAKE_GENERATOR)" -DWANT_MKMERGERMOD=0 -DWANT_LIBRARY=1

$(EVTWINLIB): $(BUILDDIR)/Makefile
	$(MAKE) -C $(BUILDDIR) -j 1 all

# Building the workers
evtwin_worker.f90: interface.py
	amusifier --type=f90 interface.py EVtwinInterface -o $@

evtwin_worker.o: evtwin_worker.f90 $(EVTWINLIB)
	$(MPIFC) -c $(FFLAGS) -I$(BUILDDIR) $<

evtwin_worker: evtwin_worker.o $(SRCDIR)/code/libint.o $(EVTWINLIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.f90
	$(FC) -c -o $@ $(FFLAGS) -I$(BUILDDIR) $<


# Which packages contain which workers?
amuse-evtwin_contains: evtwin_worker


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
	rm -rf $(BUILDDIR)
	rm -rf *.o *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

