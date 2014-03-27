AMUSE_DIR?=../..
include ${AMUSE_DIR}/config.mk

MPICC ?= mpicc

LIBNAME = libsimplehash.so

OBJS = simple_hash.o

all: test_hash test $(LIBNAME)

CFLAGS += -std=c99 -O2 -fPIC

$(LIBNAME): $(OBJS)
	$(CC) -shared -o $@ $^

test_hash: test.cpp $(OBJS)
	$(CXX) -O2 -o $@ $^

test: test_hash
	python test.py ./test_hash 123456

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
	
clean:
	rm -Rf *.o *.so test_hash
