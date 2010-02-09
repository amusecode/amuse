
CC = libtool --mode=compile gcc

all:libg6.a

libg6.a: g6lib.o
	libtool --mode=link gcc -O2 -o libg6.a g6lib.lo


clean:
	rm -Rf *.o *.lo *.a
