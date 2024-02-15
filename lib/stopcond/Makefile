HEADERS = stopcond.h stoppingconditions.mod
DYNAMIC_LIB = libstopcond.so
STATIC_LIB = libstopcond.a
PKG_CONFIG_FILE = libstopcond.pc

# Main targets
.PHONY: all
all: $(DYNAMIC_LIB) $(STATIC_LIB)

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
CFLAGS += -O2 -fPIC

OBJS = stopcond.o stoppingconditions.o

$(DYNAMIC_LIB): $(OBJS)
	$(CC) -shared -o $@ $^

$(STATIC_LIB): $(OBJS)
	$(AR) -r $(STATIC_LIB) $^
	$(RANLIB) $(STATIC_LIB)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.o %.mod &: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<


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

$(PKG_CONFIG_FILE):
	$(file >$@,$(PKG_CONFIG_CONTENTS))


# MPI support
ifdef MPICC

DYNAMIC_LIB_MPI = libstopcondmpi.so
STATIC_LIB_MPI = libstopcondmpi.a
PKG_CONFIG_FILE_MPI = libstopcondmpi.pc

all: $(DYNAMIC_LIB_MPI) $(STATIC_LIB_MPI)

.PHONY: install_mpi
install_mpi: all $(PKG_CONFIG_FILE_MPI)
	$(INSTALL) -m 644 $(STATIC_LIB_MPI) $(PREFIX)/lib/$(STATIC_LIB_MPI)
	$(INSTALL) -m 644 $(DYNAMIC_LIB_MPI) $(PREFIX)/lib/$(DYNAMIC_LIB_MPI)
	$(INSTALL) -m 644 $(PKG_CONFIG_FILE_MPI) $(PREFIX)/lib/pkgconfig/$(PKG_CONFIG_FILE_MPI)

install: install_mpi

OBJS_MPI = stopcond.mo stoppingconditions.mo

$(DYNAMIC_LIB_MPI): $(OBJS_MPI)
	$(MPICC) $(LDFLAGS) -shared -o $@ $< $(MPILIBS) $(LIBS)

$(STATIC_LIB_MPI): $(OBJS_MPI)
	$(AR) -r $(STATIC_LIB_MPI) $^
	$(RANLIB) $(STATIC_LIB_MPI)

%.mo: %.c
	$(MPICC) $(CFLAGS) -DMPILIB -c -o $@ $<

%.mo %.mod: %.f90
	$(MPIFC) $(FCFLAGS) -c -o $@ $<

define PKG_CONFIG_CONTENTS_MPI :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: stopcondmpi
Description: The stopping conditions library included with AMUSE, MPI version
Cflags: -I{includedir}
Libs: -L$${libdir} -lstopcondmpi

endef

$(PKG_CONFIG_FILE_MPI):
	$(file >$@,$(PKG_CONFIG_CONTENTS_MPI))

endif


# Clean up
.PHONY: clean
clean:
	rm -rf *.o *.mo *.mod *.a *.so *.pc
	rm -f support/config.mk support/config.log support/config.status

.PHONY: distclean
distclean: clean
	rm -rf support/autom4te.cache support/configure~

