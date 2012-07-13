MPICXX ?= mpicxx
MPICC ?= mpicc

ifneq (,$(findstring xlf, $(notdir $(FORTRAN))))
CCFLAGS += -qlanglvl=stdc99
else
CCFLAGS += -std=c99
endif
CXXFLAGS += 
LIBS += -lm
INCLUDE =

# set OpenCL and STDCL libraries
OPENCL = /home/inti/libraries/cuda
STDCL = /home/inti/libraries/coprthr_old

#---------------------------------------------

OBJS = interface.o evolve.o evolve_shared.o evolve_sf.o evolve_cc.o \
  evolve_ok.o evolve_kepler.o universal_variable_kepler.o evolve_bs.o

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all:	huayno_worker

clean:
	rm -f *.o *.bck *.pyc *.clh worker_code.cc worker_code.h worker_code-sockets.cc

distclean: clean
	rm -f huayno_worker huayno_worker_cl huayno_worker_mp huayno_worker_sockets

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py HuaynoInterface -o $@
	
worker_code-sockets.cc: interface.py
	$(CODE_GENERATOR) --type=c --mode=sockets interface.py HuaynoInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py HuaynoInterface -o $@

huayno_worker: clean __init__.py worker_code.cc worker_code.h $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) worker_code.cc $(OBJS) -o $@ $(LIBS)
	
huayno_worker_sockets: clean __init__.py worker_code-sockets.cc worker_code.h $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) worker_code-sockets.cc $(OBJS) -o $@ $(LIBS)

huayno_worker_mp: CCFLAGS   += -fopenmp
huayno_worker_mp: CXXFLAGS   += -fopenmp
huayno_worker_mp: clean __init__.py worker_code.cc worker_code.h $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) worker_code.cc $(OBJS) -o $@ $(LIBS)

huayno_worker_cl: CCFLAGS   += -fopenmp -DEVOLVE_OPENCL
huayno_worker_cl: LIBS += -L$(OPENCL)/lib/x86_64 -lOpenCL -lpthread -ldl -L$(STDCL)/lib -lstdcl
huayno_worker_cl: INCLUDE += -I$(OPENCL)/include -I$(STDCL)/include
huayno_worker_cl: clean __init__.py worker_code.cc worker_code.h evolve_kern.clh evolve_cl.o $(OBJS) 
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) worker_code.cc evolve_cl.o $(OBJS) -o $@ $(LIBS)

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $< 

.c.o: $<
	$(MPICC) $(CCFLAGS) $(INCLUDE) -c -o $@ $< 

%.clh: %.cl
	awk 'BEGIN{print "const char srcstr[]=" } {print "\""$$0"\\n\""} END{print ";"}' $< > $@
