
CFLAGS   += -Wall -g -DTOOLBOX  $(MUSE_INCLUDE_DIR)
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

BHTCDIR		= src
BHTC		= $(BHTCDIR)/BHTC
BHTC_EXTRAS	= $(BHTCDIR)/nbody.o


OBJS = bhtree_code.o 

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py


CODELIB = $(BHTCDIR)/libnearestneighbor.a

all:	compile worker_code

compile: $(BHTC)

$(BHTC):
	make -C $(BHTCDIR) all

$(CODELIB):
	make -C $(BHTCDIR) all

clean:
	rm -f *.so *.o *.pyc muse_worker.cc worker_code.cc worker_code.h 
	$(RM) muse_worker *~ worker_code *~
	make -C $(BHTCDIR) clean

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py BHTreeInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py BHTreeInterface -o $@

worker_code: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	mpicxx $(CXXFLAGS) $@.cc $(OBJS) $(CODELIB) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
	
.C.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
