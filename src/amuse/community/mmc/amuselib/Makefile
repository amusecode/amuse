MPICXX   ?= mpicxx
CXX = $(MPICXX)

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

export CFLAGS CXXFLAGS LDFLAGS

OBJS = interface.o
CODELIB = src/libsupport.a
AMUSE_DIR?=../../../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: support_worker

$(CODELIB):
	make -C src all 

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py supportInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py supportInterface -o $@

support_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(CXX) $(CXXFLAGS) $< $(OBJS) -o $@ -lm -L./src -lsupport

clean:
	rm -f *.o support_worker worker_code.cc worker_code.h
	make -C src clean	

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
