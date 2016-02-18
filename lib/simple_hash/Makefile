AMUSE_DIR?=../..
include ${AMUSE_DIR}/config.mk

MPICC ?= mpicc

LIBNAME = libsimplehash.so
STATIC_LIBNAME = libsimple_hash.a

OBJS = simple_hash.o

all: test_hash test $(LIBNAME) $(STATIC_LIBNAME)

CFLAGS += -std=c99 -O2 -fPIC

$(LIBNAME): $(OBJS)
	$(CC) -shared -o $@ $^
	
$(STATIC_LIBNAME): $(OBJS)
	ar -r  $(STATIC_LIBNAME) $^
	ranlib $(STATIC_LIBNAME)

test_hash: test.cpp $(OBJS)
	$(CXX) -O2 -o $@ $^

test: test_hash
	python test.py ./test_hash 123456

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
	
clean:
	rm -Rf *.o *.so test_hash *.a
