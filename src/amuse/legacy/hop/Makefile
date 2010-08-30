CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libhop.a

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: worker_code 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py HopInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py HopInterface -o $@

worker_code: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	mpicxx $(CXXFLAGS) $@.cc $(OBJS) $(CODELIB) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
