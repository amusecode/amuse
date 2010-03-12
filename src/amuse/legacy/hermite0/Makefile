
CFLAGS ?= -Wall -DTOOLBOX  $(MUSE_INCLUDE_DIR)
CXXFLAGS ?= $(CFLAGS) 

LDFLAGS = -lm $(MUSE_LD_FLAGS)

A_OBJS = hermite_code.o

AMUSE_DIR ?= ../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: worker_code

cleanall: clean 
	$(RM) muse_worker *~ worker_code

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py HermiteInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py HermiteInterface -o $@

worker_code:	worker_code.cc worker_code.h $(A_OBJS)
	mpicxx $(CXXFLAGS) $@.cc $(A_OBJS) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

.C.o:
	$(CC) $(CFLAGS) -c -o $@ $?

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h
