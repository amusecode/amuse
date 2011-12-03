# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

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

OBJS = libraryInterface.o

CODELIB = src/libbonsai.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: compile bonsai_worker 

$(CUDA_TK):
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "bonsai code"
	@echo "------------"
	@echo "CUDA_TK variable is not set to a valid path,"
	@echo "please set the CUDA_TK variable to the directory"
	@echo "where you installed CUDA"
	@echo "the CUDA_TK directory must contain a bin directory with the <nvcc> executable"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error

$(CUDA_SDK):
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "bonsai code"
	@echo "------------"
	@echo "CUDA_SDK variable is not set to a valid path,"
	@echo "please set the CUDA_SDK variable to the directory"
	@echo "where you installed the CUDA SDK"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error

compile: $(CUDA_TK) $(CUDA_SDK) $(CODELIB) 

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
