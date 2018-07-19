# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(strip $(AMUSE_DIR)),)
    AMUSE_DIR = ../../../..
endif
-include ${AMUSE_DIR}/config.mk

ifeq ($(strip $(CC)),)
    CC = gcc
endif
ifeq ($(strip $(CFLAGS)),)
    CFLAGS = -O3
endif

CFLAGS +=  -std=c99 -Wall -W -O3

all: mameclot_worker

mameclot_worker:
	make -C src CC="$(CC)" CFLAGS="$(CFLAGS)"  LDFLAGS="$(LDFLAGS)" F77="$(FC)" FFLAGS="$(FFLAGS)" F77FLAGS="$(F77FLAGS)"
	cp src/mameclot ./mameclot_worker

clean:
	make -C src clean
	rm -f mameclot_worker
	rm -f *.pyc *.o
