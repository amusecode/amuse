MPICXX ?= mpicxx
MPICC ?= mpicc

CXXFLAGS   += -O2 -Wall
LIBS += -lm
INCLUDE =

# set OpenCL and STDCL libraries
OPENCL = /home/inti/troep/ati-stream-sdk-v2.3-lnx64
STDCL = /home/inti/libraries/coprthr

#---------------------------------------------

OBJS = interface.o evolve.o 

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all:	worker

clean:
	rm -f *.o *.bck *.pyc *.clh worker.cc worker.h

distclean: clean
	rm -f worker worker_cl worker_mp __init__.py

worker.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py HuaynoInterface -o $@

worker.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py HuaynoInterface -o $@

worker: __init__.py worker.cc worker.h $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) worker.cc $(OBJS) -o $@ $(LIBS)

worker_mp: CXXFLAGS   += -fopenmp
worker_mp: __init__.py worker.cc worker.h $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) worker.cc $(OBJS) -o $@ $(LIBS)

worker_cl: CXXFLAGS   += -fopenmp -DEVOLVE_OPENCL
worker_cl: LIBS += -L$(OPENCL)/lib/x86_64 -lOpenCL -lpthread -ldl -L$(STDCL)/lib -lstdcl
worker_cl: INCLUDE += -I$(OPENCL)/include -I$(STDCL)/include
worker_cl: __init__.py worker.cc worker.h evolve_kern.clh evolve_cl.o $(OBJS) 
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) worker.cc evolve_cl.o $(OBJS) -o $@ $(LIBS)

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $< 

.c.o: $<
	$(MPICC) $(CXXFLAGS) $(INCLUDE) -c -o $@ $< 

%.clh: %.cl
	awk 'BEGIN{print "const char srcstr[]=" } {print "\""$$0"\\n\""} END{print ";"}' $< > $@

__init__.py:
	touch $@
