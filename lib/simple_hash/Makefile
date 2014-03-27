AMUSE_DIR?=../..
include ${AMUSE_DIR}/config.mk

MPICC ?= mpicc

LIBNAME = libsimplehash.so

OBJS = simple_hash.o

all:$(LIBNAME)

CFLAGS += -std=c99 -O2 -fPIC

$(LIBNAME): $(OBJS)
	$(CC) -shared -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
	
clean:
	rm -Rf *.o *.so
