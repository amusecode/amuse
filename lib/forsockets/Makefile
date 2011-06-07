LIBNAME = libforsockets.a
OBJS = forsockets.o

all:$(LIBNAME)

CFLAGS += -O2

$(LIBNAME): $(OBJS)
	ar -r  $(LIBNAME) $(OBJS)
	ranlib $(LIBNAME)

%.mo: %.c
	$(MPICC) $(CFLAGS) -c -o $@ $<

clean:
	rm -Rf *.o *.lo *.a *.mo
