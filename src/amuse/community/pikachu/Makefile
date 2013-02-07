# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX   ?= mpicxx
CXXFLAGS ?= -Wall -g -O2

MAP_SEQ = -DTREE_MAP
export MAP_SEQ
CXXFLAGS += $(MAP_SEQ)
SEQ_LIB ?= src/sequoia
SEQUOIA_LIBS ?= -L$(SEQ_LIB) -lsequoia
SEQUOIA_INC = -DSEQUOIA -I$(SEQ_LIB)/include

CUDA_LIBS ?= -lcuda -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64 -lcudart
CUDA_INC = -I$(CUDA_TK)/include

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CODELIB_GPU = src/libpikachu_gpu.a
CODELIB_GPU_LARGE_N = src/libpikachu_gpu_large_n.a

OBJ = interface.o
LARGE_N_OBJ = interface_large_n.o

all: $(CUDA_TK) pikachu_worker 

$(CUDA_TK):
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "pikachu code"
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
	rm -f *.o *.pyc *.ptx worker_code.cc worker_code.h 
	rm -f pikachu_worker pikachu_worker_large_n
	make -C $(SEQ_LIB) clean
	make -C src clean

$(CODELIB_GPU):
	make -C $(SEQ_LIB) all
	make -C src all
$(CODELIB_GPU_LARGE_N): $(CODELIB_GPU)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py PikachuInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py PikachuInterface -o $@

pikachu_worker: worker_code.cc worker_code.h $(CODELIB_GPU) $(OBJ)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) -I./src $< -o $@ $(OBJ) $(CODELIB_GPU) $(SC_CLIBS) $(SEQUOIA_LIBS) $(CUDA_LIBS)

pikachu_worker_large_n: worker_code.cc worker_code.h $(CODELIB_GPU_LARGE_N) $(LARGE_N_OBJ)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) -I./src $< -o $@ $(LARGE_N_OBJ) $(CODELIB_GPU_LARGE_N) $(SC_CLIBS) $(SEQUOIA_LIBS) $(CUDA_LIBS)

interface.o: interface.cc interface.h
	$(MPICXX) $(CXXFLAGS) -DSMALL $(SEQUOIA_INC) $(CUDA_INC) $(SC_FLAGS) -I./src -c -o $@ $< 

interface_large_n.o: interface.cc interface.h
	$(MPICXX) $(CXXFLAGS) $(SEQUOIA_INC) $(CUDA_INC) $(SC_FLAGS) -I./src -c -o $@ $< 
