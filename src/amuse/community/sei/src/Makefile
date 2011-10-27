all: libsei.a

libsei.a: sei.o
	rm -f libsei.a
	ar -r libsei.a sei.o
	ranlib libsei.a

clean: 
	rm libsei.a
	rm sei.o

sei.o: sei.c sei.h
	$(CC) $(CFLAGS) sei.c -c -o sei.o -lm
