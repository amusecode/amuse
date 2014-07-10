# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../..
-include ${AMUSE_DIR}/config.mk


CXX ?= g++
CC  ?= gcc
LD  ?= g++ 

.SUFFIXES: .o .cpp .ptx .cu

CUDA_TK  ?= /usr/local/cuda

OFLAGS = -g -Wall -fopenmp 

CXXFLAGS =  -fPIC $(OFLAGS) -I$(CUDA_TK)/include  -D__INCLUDE_KERNELS__ $(CL_CFLAGS)

NVCC      = $(CUDA_TK)/bin/nvcc  
NVCCFLAGS ?= -arch=sm_20

export NVCCFLAGS NVCC CXXFLAGS

# Use with Mac OS X
# NVCCFLAGS = -arch sm_12 -Xcompiler="-Duint=unsigned\ int"

CODEDIR=src/sapporo2-master

ifeq (yes,$(CUDA_ENABLED))
all:$(CODEDIR)/lib/Makefile
	make -C $(CODEDIR)/lib
	cp $(CODEDIR)/lib/libsapporo.a .
	
$(CODEDIR)/lib/Makefile:download.py
	$(PYTHON) download.py
endif

clean:
	echo $(CUDA_ENABLED)
	-rm -f libsapporo.a
	-make -C $(CODEDIR)/lib clean
	
distclean:
	-rm -f libsapporo.a
	rm -Rf src
	rm -Rf src.*
