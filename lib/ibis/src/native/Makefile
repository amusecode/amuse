MPICXX ?= mpicxx

CFLAGS   += -Wall -g -DTOOLBOX  $(MUSE_INCLUDE_DIR)
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

BHTCDIR		= src
BHTC		= $(BHTCDIR)/BHTC
BHTC_EXTRAS	= $(BHTCDIR)/nbody.o

OBJS = bhtree_code.o 

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcond
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

IBISINC ?= -I$(AMUSE_DIR)/lib/ibis/include
JNIINC ?= -I/System/Library/Frameworks/JavaVM.framework/Headers

CODELIB = $(BHTCDIR)/libbhtree.a

all:	compile bhtree_worker

compile: $(BHTC)

$(BHTC):
	make -C $(BHTCDIR) all

$(CODELIB):
	make -C $(BHTCDIR) all

clean:
	rm -f *.so *.o *.pyc muse_worker.cc worker_code.cc worker_code.h *.jnilib
	$(RM) bhtree_worker *~ worker_code *~
	make -C $(BHTCDIR) clean

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py BHTreeInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H -i amuse.support.codes.stopping_conditions.StoppingConditionInterface  interface.py BHTreeInterface -o $@

bhtree_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(SCINC) worker_code.cc $(OBJS) $(CODELIB) -o $@ $(SCLIBS)

worker_code_ibis.cc.o: worker_code_ibis.cc worker_code.h interface.py
	$(CXX) $(CXXFLAGS) $(SCINC) $(IBISINC) $(JNIINC) -c -o $@ worker_code_ibis.cc

libbhtree_worker.jnilib: worker_code_ibis.cc.o $(OBJS)
	$(CXX) -bundle -framework JavaVM -o $@ worker_code_ibis.cc.o $(OBJS) $(CODELIB) $(SCLIBS)

.cc.o: $<
	$(CXX) $(CXXFLAGS) $(SCINC) -c -o $@ $< 

.C.o: $<
	$(CXX) $(CXXFLAGS) $(SCINC) -c -o $@ $<
