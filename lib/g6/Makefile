HEADERS = g6lib.h
DYNAMIC_LIB = libg6.so
STATIC_LIB = libg6.a
PKG_CONFIG_FILE = libg6.pc


# Main targets
.PHONY: all
all: $(HEADERS) $(DYNAMIC_LIB) $(STATIC_LIB)

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
CFLAGS += -std=c99 -O2 -fPIC

OBJS = g6lib.o

$(DYNAMIC_LIB): $(OBJS)
	$(CC) -shared -o $@ $^

$(STATIC_LIB): $(OBJS)
	$(AR) -r $(STATIC_LIB) $^
	$(RANLIB) $(STATIC_LIB)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<


# Define the pkg-config file
define PKG_CONFIG_CONTENTS :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: g6lib
Description: GRAPE6 emulation library
Cflags: -I{includedir}
Libs: -L$${libdir} -lg6

endef

$(PKG_CONFIG_FILE):
	$(file >$@,$(PKG_CONFIG_CONTENTS))


# Build and run tests
test_hash: test_hash.o $(STATIC_LIB)
	$(CC) -o $@ $^

test: test_hash
	$(PYTHON) test.py ./test_hash 123456


# Clean up
.PHONY: clean
clean:
	rm -rf *.o *.a *.so *.pc test_hash
	rm -f support/config.mk support/config.log support/config.status

.PHONY: distclean
distclean: clean
	rm -rf support/autom4te.cache support/configure~

