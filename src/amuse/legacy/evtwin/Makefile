
ifndef FORTRAN
	FORTRAN = gfortran
endif

LIBDIR=src/lib
OBJDIR=src/obj

LIBTWINA=$(LIBDIR)/libtwin.a
OBJ=$(OBJDIR)/twinmuse.o 
OBJ_EVT=obj/twinmuse.o 

EVTWINLIB=src/lib/libtwin.a

all:    muse_worker

clean:
	make -C src clean
	/bin/rm -f *~ *.so *.o src/*.o src/sse muse_worker muse_worker.cc

distclean:
	make -C src distclean
	/bin/rm -f *.pyc
	
muse_worker: muse_worker.cc muse_worker.h
	make -C src all FORT=$(FORTRAN)
	make -C src $(OBJ_EVT) FORT=$(FORTRAN)
	mpicxx $< $(EVTWINLIB) $(OBJ) -o $@ -lgfortran


muse_worker.h: muse_stellar_mpi.py
	../../../../bin/create_c_header.py muse_stellar_mpi.py EVtwin > $@

muse_worker.cc: muse_stellar_mpi.py
	../../../../bin/create_c_worker.py muse_stellar_mpi.py EVtwin > $@
	
.f.o: $<
	$(FORTRAN) -c -o $@ $< 
