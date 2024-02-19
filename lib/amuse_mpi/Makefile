.PHONY: Makefile support/shared/lib-configure.mk
include support/shared/lib-configure.mk


HEADERS_MPI = amuse_mpi.h
STATIC_LIB_MPI = libamuse_mpi.a
PKG_CONFIG_FILE_MPI = libamuse_mpi.pc

CFLAGS += -g -O2 -std=c99

OBJS = amuse_mpi.mo

define PKG_CONFIG_CONTENTS_MPI :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: amuse_mpi
Description: The amuse_mpi library
Cflags: -I{includedir}
Libs: -L$${libdir} -lamuse_mpi

endef


.PHONY: support/shared/lib-targets.mk
include support/shared/lib-targets.mk

