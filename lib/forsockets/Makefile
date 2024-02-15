HEADERS = forsockets.h fortransocketsinterface.mod
DYNAMIC_LIB = libforsockets.so
STATIC_LIB = libforsockets.a
PKG_CONFIG_FILE = libforsockets.pc

# Main targets
.PHONY: all
all: $(DYNAMIC_LIB) $(STATIC_LIB)

.PHONY: install
install: all $(PKG_CONFIG_FILE)
	$(INSTALL) -m 644 $(HEADERS) $(PREFIX)/include/
	$(INSTALL) -m 644 $(STATIC_LIB) $(PREFIX)/lib/$(STATIC_LIB)
	$(INSTALL) -m 644 $(DYNAMIC_LIB) $(PREFIX)/lib/$(DYNAMIC_LIB)
	$(INSTALL) -D -m 644 $(PKG_CONFIG_FILE) $(PREFIX)/lib/pkgconfig/$(PKG_CONFIG_FILE)


# Run configure if needed
ifeq (,$(filter clean distclean, $(MAKECMDGOALS)))

include support/config.mk

support/config.mk: support/config.mk.in
	cd support && ./configure

endif


# Build the library
CFLAGS += -O2 -fPIC

OBJS = forsockets.o fortransocketsinterface.o

$(DYNAMIC_LIB): $(OBJS)
	$(CC) -shared -o $@ $^

$(STATIC_LIB): $(OBJS)
	$(AR) -r $(STATIC_LIB) $^
	$(RANLIB) $(STATIC_LIB)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.o %.mod &: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<


# Define the pkg-config file
define PKG_CONFIG_CONTENTS :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: forsockets
Description: The Fortran sockets library included with AMUSE
Cflags: -I{includedir}
Libs: -L$${libdir} -lforsockets

endef

$(PKG_CONFIG_FILE):
	$(file >$@,$(PKG_CONFIG_CONTENTS))


# Clean up
.PHONY: clean
clean:
	rm -rf *.o *.mo *.mod *.a *.so *.pc
	rm -f support/config.mk support/config.log support/config.status

.PHONY: distclean
distclean: clean
	rm -rf support/autom4te.cache support/configure~

