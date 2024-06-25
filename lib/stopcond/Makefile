.PHONY: Makefile support/shared/lib-configure.mk
include support/shared/lib-configure.mk


HEADERS = stopcond.h stopcond.inc stoppingconditions.f90
DYNAMIC_LIB = libstopcond.so
PKG_CONFIG_FILE = stopcond.pc

DYNAMIC_LIB_MPI = libstopcondmpi.so
PKG_CONFIG_FILE_MPI = stopcondmpi.pc

CFLAGS_MPI += -DMPILIB
FCFLAGS_MPI += -DMPILIB

OBJS = stopcond.o
OBJS_MPI = stopcond.mo


# Define the pkg-config file
define PKG_CONFIG_CONTENTS :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: stopcond
Description: The stopping conditions library included with AMUSE
Version: $(AMUSE_VERSION)
Cflags: -I$${includedir}
Libs: -L$${libdir} -lstopcond

endef


define PKG_CONFIG_CONTENTS_MPI :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: stopcondmpi
Description: The stopping conditions library included with AMUSE, MPI version
Version: $(AMUSE_VERSION)
Cflags: -I$${includedir}
Libs: -L$${libdir} -lstopcondmpi

endef


.PHONY: support/shared/lib-targets.mk
include support/shared/lib-targets.mk

