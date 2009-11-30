
CFLAGS = -Wall -g -DTOOLBOX  $(MUSE_INCLUDE_DIR)
CXXFLAGS = $(CFLAGS) 

LDFLAGS = -lm $(MUSE_LD_FLAGS)

BHTCDIR		= src
BHTC		= $(BHTCDIR)/BHTC
BHTC_EXTRAS	= $(BHTCDIR)/nbody.o

OBJS = muse_dynamics.o \
	$(BHTCDIR)/BHtree.o $(BHTCDIR)/gravity.o \
	$(BHTCDIR)/pgetopt.o $(BHTCDIR)/second.o $(BHTC_EXTRAS)

A_OBJS = BHTree_code.o \
	$(BHTCDIR)/BHtree.o $(BHTCDIR)/gravity.o \
	$(BHTCDIR)/pgetopt.o $(BHTCDIR)/second.o $(BHTC_EXTRAS)

all:	compile muse_worker worker_code

compile: $(BHTC)

$(BHTC):
	cd $(BHTCDIR) && $(MAKE)
	
clean:
	rm -f *.so *.o *.pyc muse_worker.cc worker_code.cc worker_code.h 
	$(RM) muse_worker *~ worker_code *~
	cd $(BHTCDIR) && $(MAKE) clean

muse_worker.cc: muse_dynamics_mpi.py
	../../../../bin/create_c_worker.py muse_dynamics_mpi.py BHTree > $@

muse_worker:	muse_worker.cc $(OBJS)
	mpicxx $@.cc $(OBJS) -o $@

worker_code.cc: interface.py
	../../../../bin/create_c_worker.py interface.py BHTree > $@

worker_code.h: interface.py
	../../../../bin/create_c_header.py interface.py BHTree > $@

worker_code: worker_code.cc worker_code.h $(A_OBJS)
	mpicxx $@.cc $(A_OBJS) -o $@

.cc.o: $<
	g++ $(CXXFLAGS) -c -o $@ $< 
.C.o: $<
	g++ $(CXXFLAGS) -c -o $@ $<
