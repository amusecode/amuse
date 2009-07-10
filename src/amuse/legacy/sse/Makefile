ifndef FORT
    FORT = gfortran
endif


TARGET = SSE_muse_interface.so muse_worker

# See if we need to run the configure script first.

MAKESYS = $(shell if [ -e makefile.sys ]; then echo 1; else echo 0; fi)
ifeq ($(MAKESYS),1)
    include makefile.sys
else
    TARGET = configure
endif

OBJ	= src/deltat.o src/evolv1.o src/hrdiag.o src/kick.o \
	  src/mlwind.o src/mrenv.o src/ran3.o src/star.o \
	  src/zcnsts.o src/zfuncs.o

compile:	$(TARGET)

configure:
	./configure.pl
	make compile

clean:
	/bin/rm -f *~ *.so *.o makefile.sys src/*.o src/sse muse_worker f2py.log

SSE_muse_interface.so:	SSE_muse_interface.f
	make -C src sse CMPLR=$(FORT)
	@echo ""
	@echo \(sending excessive f2py log output to f2py.log...\)
	@echo ""
	f2py --quiet --fcompiler=$(FORT) -c SSE_muse_interface.f \
		-m SSE_muse_interface $(OBJ) > f2py.log

muse_worker:	muse_worker.cc  SSE_muse_interface.o  $(OBJ)
	mpicxx $@.cc  SSE_muse_interface.o $(OBJ) -o $@ -lgfortran
	
.f.o: $<
	$(FORT) -c -o $@ $< 
