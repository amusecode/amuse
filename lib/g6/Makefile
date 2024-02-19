.PHONY: Makefile support/shared/lib-configure.mk
include support/shared/lib-configure.mk


HEADERS = g6lib.h
DYNAMIC_LIB = libg6.so
STATIC_LIB = libg6.a
PKG_CONFIG_FILE = libg6.pc

CFLAGS += -std=c99 -g -O2

OBJS = g6lib.o

define PKG_CONFIG_CONTENTS :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: g6lib
Description: GRAPE6 emulation library
Cflags: -I{includedir}
Libs: -L$${libdir} -lg6

endef


.PHONY: support/shared/lib-targets.mk
include support/shared/lib-targets.mk

