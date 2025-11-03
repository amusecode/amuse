# Support logic for Conda and virtualenv environments

# Detect environment
ENV_TYPE :=

ifneq (,$(VIRTUAL_ENV))
ENV_TYPE := virtualenv
ENV_NAME := $(VIRTUAL_ENV)
ENV_LIBRARY_PATH := $(VIRTUAL_ENV)/lib

HAVE_WHEEL := $(shell pip list | grep '^wheel')
HAVE_PIP := $(shell pip list | grep '^pip')
endif

ifneq (,$(CONDA_DEFAULT_ENV))
ENV_TYPE := conda
ENV_NAME := $(CONDA_DEFAULT_ENV)
ENV_LIBRARY_PATH := $(CONDA_PREFIX)/lib)

HAVE_PYPI_WHEEL := $(shell echo "$(CONDA_LIST)" | tr '^' '\n' | grep pypi | grep '^wheel')
HAVE_PYPI_PIP := $(shell echo "$(CONDA_LIST)" | tr '^' '\n' | grep pypi | grep '^pip')

HAVE_WHEEL := $(shell echo "$(CONDA_LIST)" | tr '^' '\n' | grep -v pypi | grep '^wheel')
HAVE_PIP := $(shell echo "$(CONDA_LIST)" | tr '^' '\n' | grep -v pypi | grep '^pip')
endif

ifneq (,$(HAVE_WHEEL))
ifneq (,$(HAVE_PIP))

HAVE_PIP_WHEEL := 1

endif
endif


ifneq (,$(ENV_TYPE))
ifeq (,$(HAVE_PIP_WHEEL))

NEED_PIP_WHEEL := 1

endif
endif

