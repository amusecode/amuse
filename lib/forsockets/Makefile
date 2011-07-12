LIBNAME = libforsockets.a
OBJS = forsockets.o

MPIF90 ?= openmpif90

all:$(LIBNAME) forsockets-interface.o

CFLAGS += -O2

$(LIBNAME): $(OBJS)
	ar -r  $(LIBNAME) $(OBJS)
	ranlib $(LIBNAME)

%.mo: %.c
	$(MPICC) $(CFLAGS) -c -o $@ $<
	
%.o: %.f90
	$(MPIF90) -c -o $@ $<

clean:
	rm -Rf *.o *.lo *.a *.mo *.mod
