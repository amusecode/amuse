# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx
CXXFLAGS ?= -Wall -O3 -ffast-math -funroll-loops -fomit-frame-pointer

SAPPORO_LIBS ?= -L$(AMUSE_DIR)/lib/sapporo_light -lsapporo
CUDA_LIBS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64 -lcudart


OBJS = mpi_interface.o external_field.o merge.o
GPU_OBJS = $(OBJS) evolve_gpu.o
EMU_OBJS = $(OBJS) evolve_emu.o sapporo2_dummy.o


AR = ar ruv
RANLIB = ranlib
RM = rm -f

all: Nbody_emu.out Nbody_gpu.out Nbody_gpu_pec2.out Nbody_emu_pec2.out libmi6.a libmi6_gpu.a
emulib:  libmi6.a Nbody_emu.out
gpulib:  libmi6_gpu.a

libmi6.a: $(EMU_OBJS) # Nbody_emu.o
	$(RM) $@
	$(AR) $@ $(EMU_OBJS) # Nbody_emu.o
	$(RANLIB) $@

libmi6_gpu.a: $(GPU_OBJS) # Nbody_gpu.o
	$(RM) $@
	$(AR) $@ $(GPU_OBJS) # Nbody_gpu.o
	$(RANLIB) $@

#~Nbody_emu.o:  Nbody.cc
#~	$(MPICXX) $(CXXFLAGS) -DEMU -DNOMAIN -c -o $@ $<
#~
#~Nbody_gpu.o:  Nbody.cc
#~	$(MPICXX) $(CXXFLAGS) -DSAP -DNOMAIN -c -o $@ $<

evolve_emu.o:  evolve.cc evolve.h
	$(MPICXX) -DEMU -c $(CXXFLAGS) -o $@ $<

evolve_gpu.o:  evolve.cc evolve.h
	$(MPICXX) -DSAP -c $(CXXFLAGS) -o $@ $<

evolve_emu_pec2.o:  evolve.cc evolve.h
	$(MPICXX) -DPECEC -DEMU -c $(CXXFLAGS) -o $@ $<

evolve_gpu_pec2.o:  evolve.cc evolve.h
	$(MPICXX) -DPECEC -DSAP -c $(CXXFLAGS) -o $@ $<


Nbody_emu.out:  Nbody.cc  $(EMU_OBJS)
	$(MPICXX) $(CXXFLAGS) -DEMU -o $@ $^

Nbody_gpu.out:  Nbody.cc  $(GPU_OBJS)
	$(MPICXX) $(CXXFLAGS) -DSAP -o $@ $^ -lcuda -fopenmp $(SAPPORO_LIBS) $(CUDA_LIBS)

Nbody_emu_pec2.out:  Nbody.cc  evolve_emu_pec2.o $(OBJS) sapporo2_dummy.o
	$(MPICXX) $(CXXFLAGS) -DPECEC -DEMU -o $@ $^

Nbody_gpu_pec2.out:  Nbody.cc  evolve_gpu_pec2.o $(OBJS) 
	$(MPICXX) $(CXXFLAGS) -DPECEC -DSAP -o $@ $^ -lcuda -fopenmp $(SAPPORO_LIBS) $(CUDA_LIBS)

merge.o:  merge.cc merge.h 
	$(MPICXX) -c $(CXXFLAGS) -o $@ $<

external_field.o:  external_field.cc external_field.h
	$(MPICXX) -c $(CXXFLAGS) -o $@ $<

mpi_interface.o:  mpi_interface.cc mpi_interface.h
	$(MPICXX) -c $(CXXFLAGS) -o $@ $<

sapporo2_dummy.o:  sapporo2_dummy.cc
	$(MPICXX) -c $(CXXFLAGS) -o $@ $<

clean:
	$(RM) *.o *.s Nbody_emu.out Nbody_gpu.out Nbody_emu_pec2.out Nbody_gpu_pec2.out libmi6.a libmi6_gpu.a



