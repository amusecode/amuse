

all:libg6.a

CFLAGS += -O2

libg6.a: g6lib.o
	ar -r  libg6.a g6lib.o
	ranlib libg6.a


clean:
	rm -Rf *.o *.lo *.a
