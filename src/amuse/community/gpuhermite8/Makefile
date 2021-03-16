# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

#CFLAGS   += -Wall -g -O3 -mfma
CFLAGS   += -Wall -g -mfma
#CXXFLAGS +=  -std=c++0x -I$(AMUSE_DIR)/lib/stopcond
CXXFLAGS +=  $(CFLAGS)
#LDFLAGS  +=  -lqd $(MUSE_LD_FLAGS) -L$(AMUSE_DIR)/lib/stopcond -lstopcond
LDFLAGS  +=   -lqd
NVCCFLAGS += -g -G -Xcompiler -fopenmp

#OBJS = constants.o interface.o gpuhermite8.o worker_code.o
OBJS = interface.o worker_code.o
#OBJS = interface.o gpuhermite8.o
#CODE_GENERATOR ?= $(AMUSE_DIR)/build.py
#CODELIB = src/libgpuhermite8.a

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

#gpuhermite8_worker: worker_code.cc worker_code.h  $(OBJS)
#	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) $(SC_FLAGS) $<  $(OBJS) -o $@ $(CUDA_LIBS) $(SC_CLIBS)

gpuhermite8_worker: $(OBJS)
	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) $(SC_FLAGS) $(MPI_CXXLIBS)  $(OBJS) $(CUDA_LIBS) $(SC_CLIBS) -o $@
#	$(CXX)  $(LIBS) $(LDFLAGS) $(SC_FLAGS) $(MPI_CXXLIBS)  $(OBJS)  $(SC_CLIBS) -o $@

#gpuhermite8_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
#	$(MPICXX) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#gpuhermite8.o: gpuhermite8.cu
#	$(NVCC) $(NVCCFLAGS) -c $< -o $@

#constants.o: constants.cu
#	$(NVCC) $(NVCCFLAGS) -c $< -o $@

interface.o: interface.cc
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) -c -o $@ $<

worker_code.o: worker_code.cc worker_code.h
	$(MPICXX) $(CXXFLAGS)  $(SC_FLAGS) -c -o $@ $<
