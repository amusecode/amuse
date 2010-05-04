CXX = mpicxx -g

LIBNAME = libstopcond.a
OBJS = stopcond.o

all:$(LIBNAME)

CFLAGS += -O2

$(LIBNAME): $(OBJS)
	ar -r  $(LIBNAME) $(OBJS)
	ranlib $(LIBNAME)


clean:
	rm -Rf *.o *.lo *.a
