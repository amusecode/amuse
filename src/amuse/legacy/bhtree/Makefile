
CFLAGS   += -Wall -g -DTOOLBOX  $(MUSE_INCLUDE_DIR)
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

BHTCDIR		= src
BHTC		= $(BHTCDIR)/BHTC
BHTC_EXTRAS	= $(BHTCDIR)/nbody.o

OBJS = muse_dynamics.o \
	$(BHTCDIR)/BHtree.o $(BHTCDIR)/gravity.o \
	$(BHTCDIR)/pgetopt.o $(BHTCDIR)/second.o $(BHTC_EXTRAS)

A_OBJS = bhtree_code.o \
	$(BHTCDIR)/BHtree.o $(BHTCDIR)/gravity.o \
	$(BHTCDIR)/pgetopt.o $(BHTCDIR)/second.o $(BHTC_EXTRAS)

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all:	compile worker_code

compile: $(BHTC)

$(BHTC):
	cd $(BHTCDIR) && $(MAKE)
	
clean:
	rm -f *.so *.o *.pyc muse_worker.cc worker_code.cc worker_code.h 
	$(RM) muse_worker *~ worker_code *~
	cd $(BHTCDIR) && $(MAKE) clean

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py BHTreeInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py BHTreeInterface -o $@

worker_code: worker_code.cc worker_code.h $(A_OBJS)
	mpicxx $(CXXFLAGS) $@.cc $(A_OBJS) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
.C.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
