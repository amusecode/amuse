
CFLAGS = -Wall -g -DTOOLBOX  $(MUSE_INCLUDE_DIR)
CXXFLAGS = $(CFLAGS) 

LDFLAGS = -lm $(MUSE_LD_FLAGS)

OBJS = muse_dynamics.o

all: muse_worker

cleanall: clean
	$(RM) muse_worker *~

muse_worker.cc: muse_dynamics_mpi.py
	../../../../bin/create_c_worker.py muse_dynamics_mpi.py Hermite > $@

muse_worker:	muse_worker.cc $(OBJS)
	mpicxx $@.cc $(OBJS) -o $@

.cc.o: $<
	g++ $(CXXFLAGS) -c -o $@ $< 
.C.o:
	$(CC) $(CFLAGS) -c -o $@ $?
	
clean:
	rm -f *.so *.o *.pyc muse_worker.cc
