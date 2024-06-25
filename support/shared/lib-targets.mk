# Main targets
.PHONY: all
all: $(DYNAMIC_LIB) $(DYNAMIC_LIB_MPI)

.DEFAULT_GOAL := all


# Note that the pkg-config files get built at install time, because PREFIX
# ends up inside of them, and it may not be defined when make all is called.
.PHONY: install
install: all $(PKG_CONFIG_FILE) $(PKG_CONFIG_FILE_MPI)
ifneq (,$(HEADERS)$(HEADERS_MPI))
	mkdir -p $(PREFIX)/include
endif
ifneq (,$(DYNAMIC_LIB)$(DYNAMIC_LIB_MPI))
	mkdir -p $(PREFIX)/lib
endif
ifneq (,$(PKG_CONFIG_FILE)$(PKG_CONFIG_FILE_MPI))
	mkdir -p $(PREFIX)/lib/pkgconfig
endif
ifdef HEADERS
	$(INSTALL) -m 644 $(HEADERS) $(PREFIX)/include/
endif
ifdef DYNAMIC_LIB
	$(INSTALL) -m 644 $(DYNAMIC_LIB) $(PREFIX)/lib/$(DYNAMIC_LIB)
endif
ifdef PKG_CONFIG_FILE
	$(INSTALL) -m 644 $(PKG_CONFIG_FILE) $(PREFIX)/lib/pkgconfig/$(PKG_CONFIG_FILE)
endif
ifdef HEADERS_MPI
	$(INSTALL) -m 644 $(HEADERS_MPI) $(PREFIX)/include/
endif
ifdef DYNAMIC_LIB_MPI
	$(INSTALL) -m 644 $(DYNAMIC_LIB_MPI) $(PREFIX)/lib/$(DYNAMIC_LIB_MPI)
endif
ifdef PKG_CONFIG_FILE_MPI
	$(INSTALL) -m 644 $(PKG_CONFIG_FILE_MPI) $(PREFIX)/lib/pkgconfig/$(PKG_CONFIG_FILE_MPI)
endif


# Build the library

CFLAGS += -fPIC
FCFLAGS += -fPIC

$(DYNAMIC_LIB): $(OBJS)
	$(CC) -shared -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(PKG_CONFIG_FILE):
	$(file >$@,$(PKG_CONFIG_CONTENTS))


ifdef DYNAMIC_LIB_MPI
$(DYNAMIC_LIB_MPI): $(OBJS_MPI)
	$(MPICC) $(LDFLAGS) -shared -o $@ $^ $(MPILIBS) $(LIBS)
endif

%.mo: %.c
	$(MPICC) $(CFLAGS) $(CFLAGS_MPI) -c -o $@ $<

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

