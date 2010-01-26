
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
	/bin/rm -f src/makefile.dep src/makefile.sys

distclean:
	make -C src distclean
	/bin/rm -f *.pyc

worker_code: mpi_code.f90
	make -C src lib/libtwin.a FORT=$(FORTRAN)
	cp $(MODDIR)/twin_library_v2.mod .
	cp $(MODDIR)/mesh.mod .
	mpif90 $< -o $@ $(EVTWINLIB)

mpi_code.f90: interface.py
	../../../../bin/create_f90_worker.py interface.py EVtwin > $@
	
.f.o: $<
	$(FORTRAN) -c -o $@ $< 
