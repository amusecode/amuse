.PHONY: Makefile support/shared/lib-configure.mk
include support/shared/lib-configure.mk


HEADERS = stopcond.h stoppingconditions.mod
DYNAMIC_LIB = libstopcond.so
STATIC_LIB = libstopcond.a
PKG_CONFIG_FILE = libstopcond.pc

DYNAMIC_LIB_MPI = libstopcondmpi.so
STATIC_LIB_MPI = libstopcondmpi.a
PKG_CONFIG_FILE_MPI = libstopcondmpi.pc

CFLAGS_MPI += -DMPILIB
FCFLAGS_MPI += -DMPILIB

OBJS = stopcond.o stoppingconditions.o
OBJS_MPI = stopcond.mo stoppingconditions.mo


# Define the pkg-config file
define PKG_CONFIG_CONTENTS :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: stopcond
Description: The stopping conditions library included with AMUSE
Cflags: -I{includedir}
Libs: -L$${libdir} -lstopcond

endef


define PKG_CONFIG_CONTENTS_MPI :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: stopcondmpi
Description: The stopping conditions library included with AMUSE, MPI version
Cflags: -I{includedir}
Libs: -L$${libdir} -lstopcondmpi

endef


.PHONY: support/shared/lib-targets.mk
include support/shared/lib-targets.mk

