
ifndef FORTRAN
	FORTRAN = gfortran
endif

LIBDIR=src/lib
OBJDIR=src/obj

LIBTWINA=$(LIBDIR)/libtwin.a
OBJ=$(OBJDIR)/twinmuse.o 
OBJ_EVT=obj/twinmuse.o 
MODDIR=src/modules
EVTWINLIB=src/lib/libtwin.a
 
FLIB= -L/software/local/intel64/compiler91/lib -lifcore -lifport

all:    src/makefile.dep src/makefile.sys worker_code

src/makefile.dep: src/Makefile
	@echo building makefile.deb
	@cd src && perl depend.pl code/*.f > makefile.dep
	@echo done
	

src/makefile.sys: src/Makefile
	@echo running confiure
	@cd src && perl configure.pl
	@echo done	
	
	
clean:
	make -C src clean
	/bin/rm -f *~ *.so *.o src/*.o src/sse 
	/bin/rm -f muse_worker muse_worker.f90 twin_library.mod
	/bin/rm -f muse_worker.cc  muse_worker_c muse_worker.h
	/bin/rm -f mpi_code.f90 worker_code twin_library_v2.mod

distclean:
	make -C src distclean
	/bin/rm -f *.pyc
	
muse_worker_c: muse_worker_c.cc muse_worker.h
	make -C src lib/libtwin.a $(OBJ_EVT) FORT=$(FORTRAN)
	mpicxx $< $(OBJ) -o $@ $(FLIB) $(EVTWINLIB)
	
muse_worker: muse_worker.f90	
	make -C src lib/libtwin.a FORT=$(FORTRAN)
	cp $(MODDIR)/twin_library.mod .
	cp $(MODDIR)/mesh.mod .
	mpif90 $< -o $@ $(EVTWINLIB)

worker_code: mpi_code.f90
	make -C src lib/libtwin.a FORT=$(FORTRAN)
	cp $(MODDIR)/twin_library_v2.mod .
	cp $(MODDIR)/mesh.mod .
	mpif90 $< -o $@ $(EVTWINLIB)
	

muse_worker.h: muse_stellar_mpi.py
	../../../../bin/create_c_header.py muse_stellar_mpi.py EVtwinC > $@

muse_worker_c.cc: muse_stellar_mpi.py
	../../../../bin/create_c_worker.py muse_stellar_mpi.py EVtwinC > $@

muse_worker.f90: muse_stellar_mpi.py
	../../../../bin/create_f90_worker.py muse_stellar_mpi.py EVtwin > $@
	
mpi_code.f90: interface.py
	../../../../bin/create_f90_worker.py interface.py EVtwin > $@
	
.f.o: $<
	$(FORTRAN) -c -o $@ $< 
