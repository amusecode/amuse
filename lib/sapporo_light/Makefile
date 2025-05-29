.PHONY: Makefile support/shared/lib-configure.mk
include support/shared/lib-configure.mk


HEADERS = sapporo.h sapporo_defs.h sapporo_multi.h g6lib.h
DYNAMIC_LIB = libsapporo.so
PKG_CONFIG_FILE = sapporo_light.pc

NVCCFLAGS += -D_CONSOLE -D_DEBUG -maxrregcount=32
CXXFLAGS += -DNGB -Wno-narrowing
LIBS += -lcudart

OBJS = sapporo.o send_fetch_data.o sapporoG6lib.o host_evaluate_gravity.o

define PKG_CONFIG_CONTENTS :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: sapporo_light
Description: The sapporo_light library
Version: $(VERSION)
Cflags: -I$${includedir}
Libs: -L$${libdir} -lsapporo

endef

# Main targets
.PHONY: all
all: $(DYNAMIC_LIB) $(DYNAMIC_LIB_MPI)

.DEFAULT_GOAL := all


# Install into an active venv or Conda env if no location is specified
ifdef VIRTUAL_ENV
PREFIX ?= $(VIRTUAL_ENV)
endif

ifdef CONDA_PREFIX
PREFIX ?= $(CONDA_PREFIX)
endif

ifndef PREFIX
ifeq ($(MAKECMDGOALS),install)
$(error PREFIX is not set and no virtualenv or Conda env is active.)
endif
endif


# Note that the pkg-config files get built at install time, because PREFIX
# ends up inside of them, and it may not be defined when make all is called.
.PHONY: install
install: all $(PKG_CONFIG_FILE)
	mkdir -p $(PREFIX)/include $(PREFIX)/lib/pkgconfig
	$(INSTALL) -m 644 $(HEADERS) $(PREFIX)/include/
	$(INSTALL) -m 644 $(DYNAMIC_LIB) $(PREFIX)/lib/$(DYNAMIC_LIB)
	$(INSTALL) -m 644 $(PKG_CONFIG_FILE) $(PREFIX)/lib/pkgconfig/$(PKG_CONFIG_FILE)


INSTALLED_HEADERS=$(foreach h,$(HEADERS),$(PREFIX)/include/$(h))

.PHONY: uninstall
uninstall:
	rm -f $(INSTALLED_HEADERS)
	rm -f $(PREFIX)/lib/$(DYNAMIC_LIB)
	rm -f $(PREFIX)/lib/pkgconfig/$(PKG_CONFIG_FILE)


# Build the library

CFLAGS += -fPIC
CXXFLAGS += -fPIC
NVCCFLAGS += --compiler-options '-fPIC'

$(DYNAMIC_LIB): $(OBJS) $(CUOBJS)
	$(CXX) -shared -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(PKG_CONFIG_FILE):
	$(file >$@,$(PKG_CONFIG_CONTENTS))


.PHONY: clean
clean:
	rm -rf *.o *.a *.so *.pc
	rm -f support/config.mk support/config.log support/config.status

.PHONY: distclean
distclean: clean
	rm -rf support/autom4te.cache support/configure~

