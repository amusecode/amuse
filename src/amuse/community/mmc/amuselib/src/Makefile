all: support.o libsupport.a

libsupport.a:	
	ar -r libsupport.a support.o
	ranlib libsupport.a

support.o: support.c
	$(CC) $(CFLAGS) support.c -c -o support.o

clean: 
	rm *.o *.a
