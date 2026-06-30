.PHONY: Makefile support/shared/lib-configure.mk
include support/shared/lib-configure.mk


# Main targets
.PHONY: all
all: sapporo2

.DEFAULT_GOAL := all

VERSION := 30ec27099b63057af5c542805fa67bb0e47267f0


# Install into an active venv or Conda env if no location is specified
ifdef VIRTUAL_ENV
PREFIX ?= $(VIRTUAL_ENV)
endif

ifdef CONDA_PREFIX
PREFIX ?= $(CONDA_PREFIX)
endif

ifndef PREFIX
ifeq ($(filter-out %install, $(MAKECMDGOALS)),)
$(error PREFIX is not set and no virtualenv or Conda env is active)
endif
endif


SRCDIR := src/sapporo2-$(VERSION)

sapporo2-$(VERSION).tar.gz:
	$(DOWNLOAD) https://github.com/treecode/sapporo2/archive/$(VERSION).tar.gz >$@

src:
	mkdir -p $@

$(SRCDIR): sapporo2-$(VERSION).tar.gz | src
	$(TAR) xzf $<
	mv sapporo2-$(VERSION) src/


.PHONY: sapporo2
sapporo2: | $(SRCDIR)
	$(MAKE) -C $(SRCDIR) all PREFIX=$(PREFIX)


PKG_CONFIG_FILE := sapporo2.pc

define PKG_CONFIG_CONTENTS :=
prefix=$(PREFIX)
includedir=$${prefix}/include
libdir=$${prefix}/lib

Name: sapporo2
Description: The sapporo2 library
Version: $(VERSION)
Cflags: -I$${includedir}
Libs: -L$${libdir} -lsapporo2

endef

$(PKG_CONFIG_FILE):
	# the file function is not available on old macOS make, so we make do with this
	printf '' >$@
	$(foreach line,$(PKG_CONFIG_CONTENTS),printf '%s\n' '$(line)' >>$@;)



# Note that the pkg-config files get built at install time, because PREFIX
# ends up inside of them, and it may not be defined when make all is called.
.PHONY: install
install: all $(PKG_CONFIG_FILE)
	$(MAKE) -C $(SRCDIR) install PREFIX=$(PREFIX)
	install -m 644 $(PKG_CONFIG_FILE) $(PREFIX)/lib/pkgconfig


.PHONY: uninstall
uninstall:
	$(MAKE) -C $(SRCDIR) uninstall PREFIX=$(PREFIX)
	rm -f $(PREFIX)/lib/pkgconfig/$(PKG_CONFIG_FILE)


.PHONY: clean
clean:
	-$(MAKE) -C $(SRCDIR) clean
	rm -f support/config.mk support/config.log support/config.status

.PHONY: distclean
distclean: clean
	rm -rf support/autom4te.cache support/configure~

