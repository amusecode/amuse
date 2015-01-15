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

ifdef DOWNLOAD_CODES
$(CODEDIR)/lib/Makefile:download.py
	$(PYTHON) download.py
else
$(CODEDIR)/lib/Makefile:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. Sapporo 2 will not be downloaded and build."
	@echo "If you do want Sapporo 2, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
endif
else
all:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "CUDA support is not enabled, sapporo 2 will not be build"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif

clean:
	echo $(CUDA_ENABLED)
	-rm -f libsapporo.a
	-make -C $(CODEDIR)/lib clean
	
distclean:
	-rm -f libsapporo.a
	rm -Rf src
	rm -Rf src.*
