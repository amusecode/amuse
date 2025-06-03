ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
VERSION = a6f4b64f852249a0b66d52aa6f276be517dccb64

seba.tar.gz:
	$(DOWNLOAD) https://github.com/amusecode/SeBa/archive/$(VERSION).tar.gz >$@

src: seba.tar.gz
	tar xf $<
	mv SeBa-$(VERSION) src
	# See https://github.com/amusecode/SeBa/pull/25
	cd src && patch -p1 <../toolbox.patch


# Building the code into a static library
CFLAGS += $(STOPCOND_CFLAGS) $(AMUSE_MPI_CFLAGS)
CXXFLAGS += -std=c++11       # Used by the native code
LDFLAGS +=
LDLIBS += $(STOPCOND_LIBS) $(AMUSE_MPI_LIBS)

CODELIBS = src/sstar/libsstar.a src/dstar/libdstar.a src/node/libnode.a src/std/libstd.a
CODELIBS += src/rdc/librdc.a


.PHONY: $(CODELIBS)
$(CODELIBS) &: | src
	cd src && $(MAKE) -j 1 lib


# Building the workers

INTERFACE_CFLAGS = $(CFLAGS) -Isrc/include -Isrc/include/star

seba_worker.h: interface.py
	amusifier --type=h interface.py SeBaInterface -o $@

seba_worker.cc: interface.py
	amusifier --type=c interface.py SeBaInterface -o $@

seba_worker.o: seba_worker.cc seba_worker.h | src
	$(MPICXX) -c -o $@ $(CFLAGS) $<

seba_worker: seba_worker.o interface.o $(CODELIBS)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)

interface.o: interface.cc
	$(CXX) -o $@ -c $(INTERFACE_CFLAGS) $<


# Which packages contain which workers?
amuse-seba_contains: seba_worker


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
	-$(MAKE) -C src clean
	rm -rf *.o *worker*

.PHONY: distclean
distclean: clean
	-$(MAKE) -C src distclean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

