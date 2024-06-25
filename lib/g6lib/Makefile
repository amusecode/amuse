.PHONY: Makefile support/shared/lib-configure.mk
include support/shared/lib-configure.mk


HEADERS = g6lib.h
DYNAMIC_LIB = libg6.so
PKG_CONFIG_FILE = g6lib.pc

CFLAGS += -std=c99
LIBS = -lm

OBJS = g6lib.o

define PKG_CONFIG_CONTENTS :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: g6lib
Description: GRAPE6 emulation library
Version: $(AMUSE_VERSION)
Cflags: -I$${includedir}
Libs: -L$${libdir} -lg6

endef


.PHONY: support/shared/lib-targets.mk
include support/shared/lib-targets.mk

