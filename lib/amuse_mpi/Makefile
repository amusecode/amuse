HEADERS = amuse_mpi.h
STATIC_LIB = libamuse_mpi.a
PKG_CONFIG_FILE = libamuse_mpi.pc


# Main targets
.PHONY: all
all: $(HEADERS) $(STATIC_LIB)

.PHONY: install
install: all $(PKG_CONFIG_FILE)
	$(INSTALL) -m 644 $(HEADERS) $(PREFIX)/include/
	$(INSTALL) -m 644 $(STATIC_LIB) $(PREFIX)/lib/$(STATIC_LIB)
	$(INSTALL) -D -m 644 $(PKG_CONFIG_FILE) $(PREFIX)/lib/pkgconfig/$(PKG_CONFIG_FILE)


# Run configure if needed
ifeq (,$(filter clean distclean, $(MAKECMDGOALS)))

include support/config.mk

support/config.mk: support/config.mk.in
	cd support && ./configure

endif


# Build the library
CFLAGS += -std=c99 -O2 -fPIC

OBJS = amuse_mpi.o

$(STATIC_LIB): $(OBJS)
	$(AR) -r $(STATIC_LIB) $^
	$(RANLIB) $(STATIC_LIB)

%.o: %.c
	$(MPICC) $(CFLAGS) -c -o $@ $<


# Define the pkg-config file
define PKG_CONFIG_CONTENTS :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: amuse_mpi
Description: The amuse_mpi library
Cflags: -I{includedir}
Libs: -L$${libdir} -lamuse_mpi

endef

$(PKG_CONFIG_FILE):
	$(file >$@,$(PKG_CONFIG_CONTENTS))


# Clean up
.PHONY: clean
clean:
	rm -rf *.o *.a *.pc
	rm -f support/config.mk support/config.log support/config.status

.PHONY: distclean
distclean: clean
	rm -rf support/autom4te.cache support/configure~

