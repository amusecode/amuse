AMUSE_DIR?=../../../..
include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx

CFLAGS   += -Wall -g -fPIC -DTOOLBOX  $(MUSE_INCLUDE_DIR)
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

BHTCDIR		= src
BHTC		= $(BHTCDIR)/BHTC
BHTC_EXTRAS	= $(BHTCDIR)/nbody.o

OBJS = bhtree_code.o 

CODE_GENERATOR = $(AMUSE_DIR)/build.py

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcond
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

IBISINC ?= -I$(AMUSE_DIR)/lib/ibis/include

CODELIB = $(BHTCDIR)/libbhtree.a

all:	compile bhtree_worker 

compile: $(BHTC)

$(BHTC):
	make -C $(BHTCDIR) all

$(CODELIB):
	make -C $(BHTCDIR) all

clean:
	rm -f *.so *.o *.pyc muse_worker.cc worker_code.cc worker_code.h *.jnilib worker_code-ibis.cc worker_code-sockets.cc
	$(RM) bhtree_worker *~ worker_code *~ worker_code-sockets.cc bhtree_worker_sockets
	make -C $(BHTCDIR) clean

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py BHTreeInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H -i amuse.support.codes.stopping_conditions.StoppingConditionInterface  interface.py BHTreeInterface -o $@

bhtree_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(SCINC) worker_code.cc $(OBJS) $(CODELIB) -o $@ $(SCLIBS)

worker_code-sockets.cc: interface.py
	$(CODE_GENERATOR) --type=c --mode=sockets interface.py BHTreeInterface -o $@

bhtree_worker_sockets: worker_code-sockets.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(SCINC) $(OBJS) $(CODELIB) -o $@ $(SCLIBS) $<

libbhtree_worker.jnilib: worker_code-ibis.cc  $(OBJS) $(CODELIB)
	$(CXX) $(CXXFLAGS) -bundle -framework JavaVM $(SCINC) $(IBISINC) $(JNI_INCLUDES) $(OBJS) $(CODELIB) $(SCLIBS)  -o $@ $<

libbhtree_worker.so: worker_code-ibis.cc worker_code.h interface.py $(OBJS) $(CODELIB) 
	$(CXX) $(CXXFLAGS) -fPIC -O3 -shared $(SCINC) $(IBISINC) $(JNI_INCLUDES) $(OBJS) $(CODELIB) $(SCLIBS)  -o $@ $<

.cc.o: $<
	$(CXX) $(CXXFLAGS) $(SCINC) -c -o $@ $< 

.C.o: $<
	$(CXX) $(CXXFLAGS) $(SCINC) -c -o $@ $<
