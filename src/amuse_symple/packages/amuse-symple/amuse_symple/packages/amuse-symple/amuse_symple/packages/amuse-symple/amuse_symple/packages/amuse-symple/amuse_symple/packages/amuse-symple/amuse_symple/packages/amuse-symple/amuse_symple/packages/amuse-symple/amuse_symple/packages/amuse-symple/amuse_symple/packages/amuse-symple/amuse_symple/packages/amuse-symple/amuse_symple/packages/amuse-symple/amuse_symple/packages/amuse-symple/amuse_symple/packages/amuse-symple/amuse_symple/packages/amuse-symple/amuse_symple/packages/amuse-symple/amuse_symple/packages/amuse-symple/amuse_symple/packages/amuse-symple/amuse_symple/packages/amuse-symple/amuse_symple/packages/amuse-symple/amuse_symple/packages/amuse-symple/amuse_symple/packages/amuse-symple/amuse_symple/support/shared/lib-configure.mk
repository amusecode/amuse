# Run configure if needed
ifeq (,$(filter clean distclean uninstall, $(MAKECMDGOALS)))

include support/config.mk

support/config.mk: support/config.mk.in
	cd support && ./configure $(CONFIGOPTS)

include support/shared/version.mk

export

# Install into an active venv or Conda env if no location is specified
ifdef VIRTUAL_ENV
PREFIX ?= $(VIRTUAL_ENV)
endif

ifdef CONDA_PREFIX
PREFIX ?= $(CONDA_PREFIX)
endif

ifndef PREFIX
ifeq ($(MAKECMDGOALS),install)
$(error PREFIX is not set and no virtualenv or Conda env is active.)
endif
endif

export PREFIX

endif

