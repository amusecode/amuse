.PHONY: Makefile support/shared/lib-configure.mk
include support/shared/lib-configure.mk


HEADERS_MPI = amuse_mpi.h
DYNAMIC_LIB_MPI = libamuse_mpi.so
PKG_CONFIG_FILE_MPI = amuse_mpi.pc

CFLAGS += -std=c99

OBJS_MPI = amuse_mpi.mo

define PKG_CONFIG_CONTENTS_MPI :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: amuse_mpi
Description: The amuse_mpi library
Version: $(AMUSE_VERSION)
Cflags: -I$${includedir}
Libs: -L$${libdir} -lamuse_mpi

endef


.PHONY: support/shared/lib-targets.mk
include support/shared/lib-targets.mk

