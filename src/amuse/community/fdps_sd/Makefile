# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../../../amuse
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libfdps_sd.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: fdps_sd_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ fdps_sd_worker worker_code.cc
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py fdps_sdInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py fdps_sdInterface -o $@

fdps_sd_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
