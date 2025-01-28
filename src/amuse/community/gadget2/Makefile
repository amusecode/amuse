ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
DEPFLAGS += $(STOPCONDMPI_CFLAGS) $(AMUSE_MPI_CFLAGS) $(GSL_FLAGS)

CFLAGS += $(DEPFLAGS)
CXXFLAGS += $(DEPFLAGS)

LDLIBS += $(STOPCONDMPI_LIBS) $(AMUSE_MPI_LIBS) $(GSL_LIBS)

CODELIB = libgadget.a

build_%:
	mkdir -p $@

build_%/makefile_options: makefile_options_% | build_%
	cp makefile_options_$* $@

build_%/$(CODELIB) build_%/allvars.o &: build_%/makefile_options FORCE
	$(MAKE) -C build_$* -f ../src/Makefile -j $(CPU_COUNT) all VPATH=../src


# Load code configuration
OPT :=
include makefile_options_normal
OPT_normal := $(OPT)

OPT :=
include makefile_options_nogravity
OPT_nogravity := $(OPT)

OPT :=
include makefile_options_periodic
OPT_periodic := $(OPT)

OPT :=
include makefile_options_periodic_nogravity
OPT_periodic_nogravity := $(OPT)


# Building the workers

gadget2_worker.h: interface.py
	amusifier --type=h interface.py Gadget2Interface -o $@

gadget2_worker.cc: interface.py gadget2_worker.h
	amusifier --type=c interface.py Gadget2Interface -o $@

gadget2_worker.o: gadget2_worker.cc gadget2_worker.h
	$(MPICXX) -c -o $@ $(CFLAGS) $<

gadget2_%_worker: gadget2_worker.o interface_%.o build_%/$(CODELIB) build_%/allvars.o
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

interface_%.o: interface.cc makefile_options_% gadget2_worker.h
	$(MPICXX) -c $(CXXFLAGS) -DTOOLBOX $(OPT_$*) $< -o $@


# Which packages contain which workers?
amuse-gadget2_contains: gadget2_normal_worker gadget2_nogravity_worker gadget2_periodic_worker gadget2_periodic_nogravity_worker


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
	rm -rf *.o *.mod *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

# .PHONY doesn't work for pattern rules, so we do it the old-fashioned way
FORCE:

