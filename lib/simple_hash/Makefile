.PHONY: Makefile support/shared/lib-configure.mk
include support/shared/lib-configure.mk

HEADERS = simple_hash.h
DYNAMIC_LIB = libsimple_hash.so
PKG_CONFIG_FILE = simple_hash.pc

CFLAGS += -std=c99

OBJS = simple_hash.o

define PKG_CONFIG_CONTENTS :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: simple_hash
Description: The simple_hash hashmap library included with AMUSE
Version: $(AMUSE_VERSION)
Cflags: -I$${includedir}
Libs: -L$${libdir} -lsimple_hash

endef


# Build and run tests
test_hash: test_hash.o $(DYNAMIC_LIB)
	$(CC) $(LDFLAGS) -L. -Wl,-rpath,$(PREFIX)/lib -o $@ $< -lsimple_hash

test: test_hash
	$(PYTHON) test.py ./test_hash 123456

.PHONY: test_clean
test_clean:
	rm -f test_hash

clean: test_clean

.PHONY: support/shared/lib-targets.mk
include support/shared/lib-targets.mk

