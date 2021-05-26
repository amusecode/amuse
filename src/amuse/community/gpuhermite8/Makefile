# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g -mfma
CXXFLAGS +=  $(CFLAGS) -fopenmp
LDFLAGS  +=   -L/usr/local/lib/ -lqd
NVCCFLAGS += -g -G -Xcompiler -fopenmp

OBJS = constants.o interface.o gpuhermite8.o worker_code.o
OBJS = interface.o worker_code.o
all: gpuhermite8_worker

clean:
	$(RM) -rf __pycache__
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) *~ gpuhermite8_worker worker_code.cc
#	make -C src clean

distclean: clean
#	make -C src distclean

$(CODELIB):
#	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py gpuhermite8Interface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py gpuhermite8Interface -o $@

gpuhermite8_worker: $(OBJS)
	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) $(SC_FLAGS) $(MPI_CXXLIBS)  $(OBJS) $(CUDA_LIBS) $(SC_CLIBS) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<

interface.o: interface.cc
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) -c -o $@ $<

worker_code.o: worker_code.cc worker_code.h
	$(MPICXX) $(CXXFLAGS)  $(SC_FLAGS) -c -o $@ $<
