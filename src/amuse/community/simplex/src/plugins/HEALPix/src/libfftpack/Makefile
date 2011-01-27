MPICXX ?= mpicxx

HEADERS= ls_fft.h

#include $(PARAMFILE)
INCDIR	= ../../include
LIBDIR	= ../../lib

INCLUDE = -I$(INCDIR)

all: libfftpack.a

%.o : %.cc
	$(MPICXX) $(CXXCFLAGS) -o $@ -c $< $(INCLUDE)

clean:
	$(RM) *.o *.a


FFTPACK_OBJ= fftpack.o bluestein.o ls_fft.o

fftpack.o: fftpack.h
bluestein.o: bluestein.h fftpack.h
ls_fft.o: ls_fft.h bluestein.h fftpack.h

libfftpack.a: $(FFTPACK_OBJ) $(HEADERS)
	$(ARCREATE) libfftpack.a $(FFTPACK_OBJ)
