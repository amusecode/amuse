
CODEPATH = src/octgrav_v1.7d
CODELIB = $(CODEPATH)/liboctgrav.a 


CUDA_SDK ?= /disks/koppoel1/CUDA23/cudasdk/C
CUDA_TK  ?= /disks/koppoel1/CUDA23/cuda

MUSE_INCLUDE_DIR = -I$(CODEPATH) -I$(CUDA_SDK)/common/inc -I$(CUDA_TK)/include 

CXXFLAGS+= $(MUSE_INCLUDE_DIR) 

OBJS = integrator.o octgrav_code.o 

AMUSE_DIR?=../../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py

CUDA_LIBDIRS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64
CUDA_LIBS ?= -lcudart

LIBS += $(CUDA_LIBDIRS) $(CUDA_LIBS) 

$(warning make does: : $(MAKE))

all:  compile worker_code

$(CUDA_TK):
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "octgrav code"
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
	@echo "octgrav code"
	@echo "------------"
	@echo "CUDA_SDK variable is not set to a valid path,"
	@echo "please set the CUDA_SDK variable to the directory"
	@echo "where you installed the CUDA SDK"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error

compile: $(CUDA_TK) $(CUDA_SDK) octgrav 

octgrav: 
	$(MAKE) -C $(CODEPATH)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py OctgravInterface -o $@ 

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py OctgravInterface -o $@

worker_code: worker_code.cc worker_code.h $(OBJS)
	mpicxx $(CXXFLAGS) $@.cc $(OBJS) $(LIBS) $(CODELIB) -o $@

clean:
	rm -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) worker_code *~
	$(MAKE) -C $(CODEPATH) clean	

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

.C.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: octgrav
