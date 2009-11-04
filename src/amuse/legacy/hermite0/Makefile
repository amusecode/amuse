
CFLAGS = -Wall -g -DTOOLBOX  $(MUSE_INCLUDE_DIR)
CXXFLAGS = $(CFLAGS) 

LDFLAGS = -lm $(MUSE_LD_FLAGS)

OBJS = muse_dynamics.o

A_OBJS = hermite_code.o

all: muse_worker worker_code

cleanall: clean 
	$(RM) muse_worker *~ worker_code

muse_worker.cc: muse_dynamics_mpi.py
	../../../../bin/create_c_worker.py muse_dynamics_mpi.py Hermite > $@
	
worker_code.cc: interface.py
	../../../../bin/create_c_worker.py interface.py Hermite > $@

worker_code.h: interface.py
	../../../../bin/create_c_header.py interface.py Hermite > $@

muse_worker:	muse_worker.cc $(OBJS)
	mpicxx $@.cc $(OBJS) -o $@

worker_code:	worker_code.cc worker_code.h $(A_OBJS)
	mpicxx $@.cc $(A_OBJS) -o $@

.cc.o: $<
	g++ $(CXXFLAGS) -c -o $@ $< 
.C.o:
	$(CC) $(CFLAGS) -c -o $@ $?
	
clean:
	rm -f *.so *.o *.pyc muse_worker.cc worker_code.cc worker_code.h
