# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX   ?= mpicxx
CXXFLAGS ?= -Wall -g -O2

MAP_SEQ ?= -DTREE_MAP



CCC = $(MPICXX) $(CXXFLAGS) -DMPICH_IGNORE_CXX_SEEK $(MAP_SEQ) $(OPENMP_CFLAGS)

CUDA_TK  ?= /usr/local/cuda
CUDA_LIBS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64 -lcudart -lcuda 

SEQ_LIBDIR = ./sequoia
LIB_SEQ = -lsequoia -L$(SEQ_LIBDIR)
INC_SEQ = -DSEQUOIA -I$(SEQ_LIBDIR)/include

#OFLAGS_SEQ = -O3 -g -Wall -fopenmp 
CFLAGS_SEQ =  -fPIC $(OFLAGS_SEQ) -I$(CUDA_TK)/include/

AR = ar ruv
RANLIB = ranlib
RM = rm -f


OBJS_SEQ = BHtree.o soft_system_seq.o hard_system_seq.o nbody_system_seq.o
OBJS_SEQ_LARGE_N = BHtree_large_n.o soft_system_seq_large_n.o hard_system_seq_large_n.o nbody_system_seq_large_n.o

CODELIB_GPU = libpikachu_gpu.a
CODELIB_GPU_LARGE_N = libpikachu_gpu_large_n.a

all: copy hybrid_seq.out $(CODELIB_GPU) $(CODELIB_GPU_LARGE_N)

copy: $(SEQ_LIBDIR)/*.ptx
	cp $(SEQ_LIBDIR)/*.ptx ../

$(CODELIB_GPU): $(OBJS_SEQ)
	$(RM) $@
	$(AR) $@ $(OBJS_SEQ)
	$(RANLIB) $@

$(CODELIB_GPU_LARGE_N): $(OBJS_SEQ_LARGE_N)
	$(RM) $@
	$(AR) $@ $(OBJS_SEQ_LARGE_N)
	$(RANLIB) $@


hybrid_seq.out: Nbody.cc $(OBJS_SEQ_LARGE_N) distribution.h const.h
	$(CCC) $(INC_SEQ) $(CFLAGS_SEQ) -o $@ $<  $(OBJS_SEQ_LARGE_N) $(LIB_SEQ) $(CUDA_LIBS)


nbody_system_seq.o: nbody_system.cc system.h const.h
	$(CCC) -DSMALL -c $(INC_SEQ) $(CFLAGS_SEQ) -o $@ $<

hard_system_seq.o: hard_system.cc system.h const.h  force.h
	$(CCC) -DSMALL -c $(INC_SEQ) $(CFLAGS_SEQ) -o $@ $<

soft_system_seq.o: soft_system.cc system.h const.h
	$(CCC) -DSMALL -c $(INC_SEQ) $(CFLAGS_SEQ) -o $@ $<

BHtree.o: BHtree.cc BHtree.h particle.h mpi_interface.h const.h
	$(CCC) -DSMALL -c -o $@ $<


nbody_system_seq_large_n.o: nbody_system.cc system.h const.h
	$(CCC) -c $(INC_SEQ) $(CFLAGS_SEQ) -o $@ $<

hard_system_seq_large_n.o: hard_system.cc system.h const.h  force.h
	$(CCC) -c $(INC_SEQ) $(CFLAGS_SEQ) -o $@ $<

soft_system_seq_large_n.o: soft_system.cc system.h const.h
	$(CCC) -c $(INC_SEQ) $(CFLAGS_SEQ) -o $@ $<

BHtree_large_n.o: BHtree.cc BHtree.h particle.h mpi_interface.h const.h
	$(CCC) -c -o $@ $<


clean:
	$(RM) *.s *.o *.out *.ptx $(CODELIB_GPU) $(CODELIB_GPU_LARGE_N)

