
# Makefile to build the worker code.  First draft was auto-generated;
# subsequent versions are hand-coded, so save a copy.  Note that
# all compiles must use MPICXX, not CXX.

# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

ALL	 = hacs64_worker
OPT	  = -I$(AMUSE_DIR)/lib/stopcond 

MPICXX   ?= mpicxx

CFLAGS   += $(OPT)
CFLAGS   += -I./src -I./src/ON_neib
CXXFLAGS += $(CFLAGS)  
LDFLAGS  += -L$(AMUSE_DIR)/lib/stopcond -lstopcond -lm $(MUSE_LD_FLAGS)

CUDA_TK  ?= /usr/local/cuda

CUDA_LIBDIRS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64
CUDA_LIBS ?= $(CUDA_LIBDIRS) -lcudart

LDFLAGS_GPU =	$(LDFLAGS) \
		$(CUDA_LIBS) \

CXXFLAGS += -I $(CUDA_TK)/include


OBJS = interface.o
SRCDIR = src
BUILDDIRCPU = build_cpu
CXXFLAGS += -I../src -I../src/ON_neib

CODELIB = src/libhacs64gpu.a
CODELIBCPU = $(BUILDDIRCPU)/libhacs64cpu.a

AMUSE_DIR ?= ../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: $(ALL)

save:
	@ for f in Makefile interface.py interface.cc src/Makefile; do \
		cp -p $$f $$f.save; done

$(CUDA_TK):
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "hacs64 code"
	@echo "------------"
	@echo "CUDA_TK variable is not set to a valid path,"
	@echo "please set the CUDA_TK variable to the directory"
	@echo "where you installed CUDA"
	@echo "the CUDA_TK directory must contain a bin directory with the <nvcc> executable"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error


$(CODELIB):
	make -C src all CXX='$(MPICXX)'


$(BUILDDIRCPU)/Makefile: $(SRCDIR)/Makefile 
	-mkdir $(BUILDDIRCPU) 
	cp $(SRCDIR)/Makefile $(BUILDDIRCPU)/Makefile

$(CODELIBCPU): $(BUILDDIRCPU)/Makefile
	make -C $(BUILDDIRCPU) cpu CXX='$(MPICXX)' VPATH=../src

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py Hacs64Interface -o $@
	
interface.h: interface.py
	$(CODE_GENERATOR) --type=H  -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py Hacs64Interface -o $@

hacs64_worker: worker_code.cc interface.h $(CODELIBCPU) $(OBJS)
	$(MPICXX) $(CXXFLAGS) worker_code.cc $(OBJS) $(CODELIBCPU) -o $@  $(LDFLAGS)  $(LIBS)

hacs64_worker_gpu:  $(CUDA_TK) worker_code.cc interface.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) worker_code.cc $(OBJS) $(CODELIB) -o $@ $(LDFLAGS_GPU) -DGPU  $(LIBS)

%.o: %.cc
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 

clean:
	$(RM) -r -f *.so *.dSYM *.o *.mpio *.pyc worker_code.cc interface.h 
	$(RM) -f *~ hacs64_worker hacs64_worker_cpu worker_code.cc
	$(RM) -f worker_code-sockets.cc
	$(RM) -rf $(BUILDDIRCPU)
	make -C src clean
