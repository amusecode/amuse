ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
VERSION = 4f97bfe11e8c638fdda744ca288e57565efe718a

tidymess.tar.gz:
	$(DOWNLOAD) https://github.com/tidymess-code/tidymess/archive/$(VERSION).tar.gz >$@

src/tidymess: tidymess.tar.gz
	tar xf $<
	mv tidymess-$(VERSION) src/tidymess

# Building the code into a static library
DEPFLAGS += $(STOPCOND_CFLAGS)
CXXFLAGS += $(DEPFLAGS)

LDLIBS += -lm $(STOPCOND_LIBS)

CODELIB = src/libtidymess.a

.PHONY: $(CODELIB)
$(CODELIB): | src/tidymess
	$(MAKE) -C src -j $(CPU_COUNT) all


# Building the workers
tidymess_worker.h: interface.py
	amusifier --type=h interface.py TidymessInterface -o $@

tidymess_worker.cc: interface.py
	amusifier --type=c interface.py TidymessInterface -o $@

tidymess_worker.o: tidymess_worker.cc tidymess_worker.h
	$(MPICXX) -c -o $@ $(CXXFLAGS) $<

tidymess_worker: interface.o tidymess_worker.o $(CODELIB)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)


interface.o: interface.cc tidymess_worker.h | src/tidymess
	$(MPICXX) -o $@ -c -I src/tidymess/integrator/include $(CXXFLAGS) $<


# Which packages contain which workers?
amuse-tidymess_contains: tidymess_worker


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
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache
