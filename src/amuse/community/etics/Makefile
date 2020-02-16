# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk
GPUARCH ?= 35
CUDAFLAGS += -arch=sm_$(GPUARCH)
######### ugly!!!! should be in the amuse base make file


MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = build/libetics.a

CODE_GENERATOR ?= $(AMUSE_DIR)/build.py

CUDA_TK ?= /usr/local/cuda/
CUDA_LIBS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64 -lcuda 

all: compile etics_worker

$(CUDA_TK):
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "Etics code"
	@echo "------------"
	@echo "CUDA_TK variable is not set to a valid path,"
	@echo "please set the CUDA_TK variable to the directory"
	@echo "where you installed CUDA"
	@echo "the CUDA_TK directory must contain a bin directory with the <nvcc> executable"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error


clean:
	$(RM) *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ etics_worker worker_code.cc
	make -C src clean

build/config.mk: 
	-mkdir build
	cd build && ../src/configure

$(CODELIB): build/config.mk
	make -C build libetics.a

compile: $(CUDA_TK) $(CODELIB)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py EticsInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py EticsInterface -o $@

etics_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) -L/usr/local/cuda/lib64 $< $(OBJS) $(CODELIB) -o $@ -lcudart
	echo the above compilation directive is not nice

.SUFFIXES: .cu .o

.cu.o: $<
	$(NVCC) -Xcompiler="$(CXXFLAGS)" -c -o $@ $<


# .cc.o: $<
# 	$(CXX) $(CXXFLAGS) -c -o $@ $<
