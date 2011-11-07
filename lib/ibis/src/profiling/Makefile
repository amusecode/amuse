MPICC ?= mpicc

LIBNAME = libibisprofiling.a
OBJS = ibisprofiling.o ibisstatistics.o

all:$(LIBNAME)

CFLAGS += -O2

$(LIBNAME): $(OBJS)
	ar -r  $(LIBNAME) $(OBJS)
	ranlib $(LIBNAME)

%.o: %.c
	$(MPICC) $(CFLAGS) -c -o $@ $<

clean:
	rm -Rf *.o *.lo *.a *.mo
