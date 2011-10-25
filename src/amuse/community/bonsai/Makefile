MPICXX   ?= openmpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CUDA_SDK ?= /disks/koppoel1/CUDA23/cudasdk/C
CUDA_TK ?= /usr/local/cuda/
CUDA_LIBDIRS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64

#CXXFLAGS += -I$(CUDA_TK)/include
LIBS = -lcuda -lOpenCL $(LDFLAGS) $(CUDA_LIBDIRS) -L/usr/lib64/nvidia
INCLUDES = -I$(CUDA_SDK)/common/inc -I$(CUDA_TK)/include 

OBJS = interface.o libraryInterface.o

CODELIB = src/libbonsai.a

#CODELIB += $(LIBS)

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: bonsai_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ bonsai_worker worker_code.cc
	make -C src clean

$(CODELIB):
	make -C src all CXX="$(MPICXX)" CUDA_TK="$(CUDA_TK)" LD="$(MPICXX)" LDFLAGS="$(LIBS)" INCLUDES="$(INCLUDES)"

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py BonsaiInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py BonsaiInterface -o $@

bonsai_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) $< $(OBJS) $(CODELIB) $(LIBS) -o $@

.cpp.o: $<
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $< 
.cc.o: $<
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $< 
