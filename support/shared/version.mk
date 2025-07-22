AMUSE_VERSION := $(patsubst v%,%,$(shell git describe --tags))

ifeq (,$(AMUSE_VERSION))
    AMUSE_VERSION := $(shell grep -v '^#' ../../VERSION)
endif

export AMUSE_VERSION

