# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

CC ?= gcc
CFLAGS ?= -O3

CFLAGS +=  -std=c99 -Wall -W -O3

mameclot_worker:
	make -C src CC="$(CC)" CFLAGS="$(CFLAGS)" F77="$(FC)" FFLAGS="$(FFLAGS)" F77FLAGS="$(F77FLAGS)"
	cp src/mameclot ./mameclot_worker

clean:
	make -C src clean
	rm -f mameclot_worker
	rm -f *.pyc *.o
