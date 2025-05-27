AMUSE_VERSION := $(patsubst v%,%,$(shell git describe --tags))

export AMUSE_VERSION

