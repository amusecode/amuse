.PHONY: Makefile support/shared/lib-configure.mk
include support/shared/lib-configure.mk


HEADERS = forsockets.h fortransocketsinterface.f90
DYNAMIC_LIB = libforsockets.so
STATIC_LIB = libforsockets.a
PKG_CONFIG_FILE = forsockets.pc

OBJS = forsockets.o

define PKG_CONFIG_CONTENTS :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: forsockets
Description: The Fortran sockets library included with AMUSE
Version: $(AMUSE_VERSION)
Cflags: -I$${includedir}
Libs: -L$${libdir} -lforsockets

endef


.PHONY: support/shared/lib-targets.mk
include support/shared/lib-targets.mk

