ifndef FORT
    FORT = gfortran
endif



OBJ	= src/deltat.o src/evolv1.o src/hrdiag.o src/kick.o \
	  src/mlwind.o src/mrenv.o src/ran3.o src/star.o \
	  src/zcnsts.o src/zfuncs.o

all:    muse_worker1


clean:
	/bin/rm -f *~ *.so *.o src/*.o src/sse muse_worker1 muse_worker.f90

muse_worker1:	muse_worker.f90  SSE_muse_interface.o  $(OBJ)
	make -C src sse CMPLR=$(FORT)
	mpif90 $^ -o $@ 


muse_worker.f90: muse_stellar_mpi.py
	../../../../bin/create_f90_worker.py muse_stellar_mpi.py SSE > $@

.f.o: $<
	$(FORT) -c -fPIC -o $@ $< 
