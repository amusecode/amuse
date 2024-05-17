# Main targets
.PHONY: all
all: $(DYNAMIC_LIB) $(STATIC_LIB) $(DYNAMIC_LIB_MPI) $(STATIC_LIB_MPI)

.DEFAULT_GOAL := all


# Note that the pkg-config files get built at install time, because PREFIX
# ends up inside of them, and it may not be defined when make all is called.
.PHONY: install
install: all $(PKG_CONFIG_FILE) $(PKG_CONFIG_FILE_MPI)
ifdef HEADERS
	$(INSTALL) -m 644 $(HEADERS) $(PREFIX)/include/
endif
ifdef STATIC_LIB
	$(INSTALL) -m 644 $(STATIC_LIB) $(PREFIX)/lib/$(STATIC_LIB)
endif
ifdef DYNAMIC_LIB
	$(INSTALL) -m 644 $(DYNAMIC_LIB) $(PREFIX)/lib/$(DYNAMIC_LIB)
endif
ifdef PKG_CONFIG_FILE
	$(INSTALL) -D -m 644 $(PKG_CONFIG_FILE) $(PREFIX)/lib/pkgconfig/$(PKG_CONFIG_FILE)
endif
ifdef HEADERS_MPI
	$(INSTALL) -m 644 $(HEADERS_MPI) $(PREFIX)/include/
endif
ifdef STATIC_LIB_MPI
	$(INSTALL) -m 644 $(STATIC_LIB_MPI) $(PREFIX)/lib/$(STATIC_LIB_MPI)
endif
ifdef DYNAMIC_LIB_MPI
	$(INSTALL) -m 644 $(DYNAMIC_LIB_MPI) $(PREFIX)/lib/$(DYNAMIC_LIB_MPI)
endif
ifdef PKG_CONFIG_FILE_MPI
	$(INSTALL) -D -m 644 $(PKG_CONFIG_FILE_MPI) $(PREFIX)/lib/pkgconfig/$(PKG_CONFIG_FILE_MPI)
endif


# Build the library

CFLAGS += -fPIC
FCFLAGS += -fPIC

$(DYNAMIC_LIB): $(OBJS)
	$(CC) -shared -o $@ $^

$(STATIC_LIB): $(OBJS)
	$(AR) -r $(STATIC_LIB) $^
	$(RANLIB) $(STATIC_LIB)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.o %.mod &: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

$(PKG_CONFIG_FILE):
	$(file >$@,$(PKG_CONFIG_CONTENTS))


ifdef DYNAMIC_LIB_MPI
$(DYNAMIC_LIB_MPI): $(OBJS_MPI)
	$(MPICC) $(LDFLAGS) -shared -o $@ $^ $(MPILIBS) $(LIBS)
endif

ifdef STATIC_LIB_MPI
$(STATIC_LIB_MPI): $(OBJS_MPI)
	$(AR) -r $(STATIC_LIB_MPI) $^
	$(RANLIB) $(STATIC_LIB_MPI)
endif

%.mo: %.c
	$(MPICC) $(CFLAGS) $(CFLAGS_MPI) -c -o $@ $<

%.mo %.mod &: %.f90
	$(MPIFC) $(FCFLAGS) $(FCFLAGS_MPI) -c -o $*.mo $<

$(PKG_CONFIG_FILE_MPI):
	$(file >$@,$(PKG_CONFIG_CONTENTS_MPI))


# Clean up
.PHONY: clean
clean:
	rm -rf *.o *.mo *.mod *.a *.so *.pc
	rm -f support/config.mk support/config.log support/config.status

.PHONY: distclean
distclean: clean
	rm -rf support/autom4te.cache support/configure~

