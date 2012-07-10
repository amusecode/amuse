# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX   ?= openmpicxx


SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcondmpi   
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CUDA_TK ?= /usr/local/cuda/
CUDA_LIBDIRS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64

#CXXFLAGS += -I$(CUDA_TK)/include
LIBS = -lcuda  $(LDFLAGS) $(CUDA_LIBDIRS) -L/usr/lib64/nvidia $(SCLIBS)
INCLUDES = -I$(CUDA_TK)/include $(SCINC)

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


compile: $(CUDA_TK) $(CODELIB) 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ bonsai_worker worker_code.cc worker_code-sockets.cc bonsai_worker_sockets
	make -C src clean

$(CODELIB):
	make -C src all CXX="$(MPICXX)" CUDA_TK="$(CUDA_TK)" LD="$(MPICXX)" LDFLAGS="$(LIBS)" INCLUDES="$(INCLUDES)"

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py BonsaiInterface -o $@

worker_code-sockets.cc: interface.py
	$(CODE_GENERATOR) --type=c --mode=sockets interface.py BonsaiInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H -i amuse.support.codes.stopping_conditions.StoppingConditionInterface  interface.py BonsaiInterface -o $@


bonsai_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) -o $@ $(SCLIBS) $(INCLUDES) $< $(OBJS) $(CODELIB) $(LIBS)

bonsai_worker_sockets: worker_code-sockets.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) -o $@ $(SCLIBS) $(INCLUDES) $< $(OBJS) $(CODELIB) $(LIBS)

.cpp.o: $<
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $< 
.cc.o: $<
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $< 
