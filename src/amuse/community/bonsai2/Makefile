# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk

MPICXX   ?= openmpicxx

CFLAGS   += -Wall -g -std=c++11 -Wno-strict-aliasing -Wno-reorder
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CUDA_TK ?= /usr/local/cuda/
CUDA_LIBS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64 -lcuda 

#CXXFLAGS += -I$(CUDA_TK)/include
LIBS += $(CUDA_LIBS) $(LDFLAGS) $(CUDA_LIBDIRS) -L/usr/lib64/nvidia $(SC_MPI_CLIBS)
INCLUDES = -I$(CUDA_TK)/include $(SC_FLAGS)

OBJS = libraryInterface.o

CODELIB = src/bonsai/runtime/libbonsai_amuse.a

DOWNLOAD_FROM_WEB = $(PYTHON) ./download.py

CODE_GENERATOR ?= $(AMUSE_DIR)/build.py

all: compile bonsai_worker 

ifdef DOWNLOAD_CODES
$(CODEDIR)/Makefile:
	make -C . download
else
$(CODEDIR)/Makefile:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. Bonsai2 will not be downloaded and built."
	@echo "If you do want Bonsai2, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif

download:
	$(RM) -Rf .pc
	$(RM) -Rf src
	mkdir src
	$(DOWNLOAD_FROM_WEB)
	$(PATCH_FILES)

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


compile: $(CUDA_TK) $(CODELIB)
	cd src/bonsai/runtime
	cmake -DUSE_MPI=0 -DUSE_MPIMT=0 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_THRUST=1 -DUSE_CUB=0
	make -j
	make

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ bonsai_worker worker_code.cc worker_code-sockets.cc bonsai_worker_sockets
	make -C src clean

#$(CODELIB):
#	make -C src all CXX="$(MPICXX)" CUDA_TK="$(CUDA_TK)" LD="$(MPICXX)" LDFLAGS="$(LIBS)" INCLUDES="$(INCLUDES)"

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py BonsaiInterface2 -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h -i amuse.support.codes.stopping_conditions.StoppingConditionInterface  interface.py BonsaiInterface2 -o $@


bonsai_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(SC_MPI_CLIBS) $(INCLUDES) $< $(OBJS) $(CODELIB) $(LIBS)

.cpp.o: $<
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $< 
.cc.o: $<
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $< 
