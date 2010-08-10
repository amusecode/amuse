LIBNAME = libstopcond.a
LIBNAMEMPI = libstopcondmpi.a
OBJS = stopcond.o
OBJSMPI = stopcond.mo

all:$(LIBNAME) $(LIBNAMEMPI)

CFLAGS += -O2

$(LIBNAME): $(OBJS)
	ar -r  $(LIBNAME) $(OBJS)
	ranlib $(LIBNAME)

$(LIBNAMEMPI): $(OBJSMPI)
	ar -r  $(LIBNAMEMPI) $(OBJSMPI)
	ranlib $(LIBNAMEMPI)


%.mo: %.c
	mpicc $(CXXFLAGS) -DMPILIB  -c -fno-leading-underscore -o $@ $< 


clean:
	rm -Rf *.o *.lo *.a
