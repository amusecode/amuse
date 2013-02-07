AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx

CFLAGS   += -DTOOLBOX  $(MUSE_INCLUDE_DIR)
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

BHTCDIR		= src
BHTC		= $(BUILDDIR)/BHTC
BHTC_EXTRAS	= $(BUILDDIR)/nbody.o

BUILDDIR = build
GPU_BUILDDIR = build_gpu

OBJS = interface.o 

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CODELIB = $(BUILDDIR)/libbhtree.a
GPUCODELIB = $(GPU_BUILDDIR)/libbhtree.a
G6LIBS ?= -L$(AMUSE_DIR)/lib/g6 -lg6

SAPPORO_LIBDIRS ?= -L$(AMUSE_DIR)/lib/sapporo_light
OPENMP_CFLAGS ?= 
SAPPORO_LIBS ?= $(SAPPORO_LIBDIRS) -lsapporo $(OPENMP_CFLAGS)
CUDA_LIBDIRS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64
CUDA_LIBS ?= $(CUDA_LIBDIRS) -lcudart
GPULIBS  =  $(CUDA_LIBS) $(SAPPORO_LIBS) -lpthread

all:	compile bhtree_worker

python_worker.o: python_worker.c
	$(CXX) $(CXXFLAGS) $(PYTHONDEV_CFLAGS) $(SC_FLAGS) -c -o $@ $< 
	
#bhtree_worker_python: python_worker.o
#	$(MPICXX) $(CXXFLAGS) $(PYTHONDEV_LIBS) $(SC_FLAGS) python_worker.o $(OBJS) $(CODELIB) -o $@ $(SC_LIBS)
    

compile: $(BHTC)

$(BUILDDIR)/Makefile: $(BHTCDIR)/Makefile 
	-mkdir $(BUILDDIR) 
	cp $(BHTCDIR)/Makefile $(BUILDDIR)/Makefile


$(GPU_BUILDDIR)/Makefile: $(BHTCDIR)/Makefile 
	-mkdir $(GPU_BUILDDIR) 
	cp $(BHTCDIR)/Makefile $(GPU_BUILDDIR)/Makefile
	
$(BHTC): $(BUILDDIR)/Makefile
	make -C $(BUILDDIR) all VPATH=../src
$(CODELIB): $(BUILDDIR)/Makefile
	make -C $(BUILDDIR) all VPATH=../src
	
$(GPUCODELIB): $(GPU_BUILDDIR)/Makefile
	make -C $(GPU_BUILDDIR) libbhtree.a CPPFLAGS=-DGPU VPATH=../src

clean:
	$(RM) -f *.so *.o *.pyc muse_worker.cc worker_code.cc worker_code.h *.jnilib worker_code-ibis.cc worker_code-sockets.cc
	$(RM) bhtree_worker *~ worker_code *~ worker_code-sockets.cc bhtree_worker_sockets
	$(RM) -rf *.dSYM
	$(RM) -Rf $(BUILDDIR) $(GPU_BUILDDIR)
	$(RM) -f bhtree_worker_g6 bhtree_worker_gpu
	make -C $(BHTCDIR) clean

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c --needs-mpi=false interface.py BHTreeInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H -i amuse.support.codes.stopping_conditions.StoppingConditionInterface  interface.py BHTreeInterface -o $@

bhtree_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) worker_code.cc $(OBJS) $(CODELIB) -o $@ $(SC_CLIBS)
	
#bhtree_worker_g6: worker_code.cc worker_code.h $(GPUCODELIB) $(OBJS)
#	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) worker_code.cc $(OBJS) $(GPUCODELIB) -o $@ $(SC_LIBS) $(G6LIBS)
	
.cc.o: $<
	$(CXX) $(CXXFLAGS) $(SC_FLAGS) -c -o $@ $< 

.C.o: $<
	$(CXX) $(CXXFLAGS) $(SC_FLAGS) -c -o $@ $<
