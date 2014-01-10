# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx
MPICC ?= mpicc

OPENMP_CFLAGS?= -fopenmp

ifneq (,$(findstring xlc, $(notdir $(CC))))
CFLAGS += -qlanglvl=extc99
else
CFLAGS += -std=gnu99
endif
CFLAGS +=
CXXFLAGS +=  
LIBS += -lm
INCLUDE =


CODE_GENERATOR = $(AMUSE_DIR)/build.py

# set OpenCL and STDCL libraries
OPENCL = /home/inti/libraries/cuda
STDCL = /home/inti/libraries/coprthr

#---------------------------------------------

OBJS = interface.o

SRCDIR = src

BUILDDIR = build
OPENCL_BUILDDIR = build_cl
OPENMP_BUILDDIR = build_mp

CL_LIBS?=-L$(OPENCL)/lib/x86_64 -lOpenCL -lpthread -ldl -L$(STDCL)/lib -lstdcl
CL_FLAGS?=-I$(OPENCL)/include -I$(STDCL)/include

all:	huayno_worker

clean:
	rm -f *.o *.bck *.pyc *.clh
	rm -Rf $(BUILDDIR)
	rm -Rf $(OPENCL_BUILDDIR)
	rm -Rf $(OPENMP_BUILDDIR)
	rm -f huayno_worker huayno_worker_cl huayno_worker_mp huayno_worker_sockets

distclean: clean

$(BUILDDIR)/Makefile: $(SRCDIR)/Makefile 
	-mkdir $(BUILDDIR) 
	cp $(SRCDIR)/Makefile $(BUILDDIR)/Makefile
	
$(BUILDDIR)/libhuayno.a: $(BUILDDIR)/Makefile $(LIBFILES) src/Makefile
	make -C $(BUILDDIR) all VPATH=../src 
	
	
$(OPENMP_BUILDDIR)/Makefile: $(SRCDIR)/Makefile 
	-mkdir $(OPENMP_BUILDDIR) 
	cp $(SRCDIR)/Makefile $(OPENMP_BUILDDIR)/Makefile
	
$(OPENMP_BUILDDIR)/libhuayno.a: $(OPENMP_BUILDDIR)/Makefile $(LIBFILES) src/Makefile
	make -C $(OPENMP_BUILDDIR) all VPATH=../src  CFLAGS="$(CFLAGS)" CXXFLAGS="$(CXXFLAGS)" 
	
$(OPENCL_BUILDDIR)/Makefile: $(SRCDIR)/Makefile 
	-mkdir $(OPENCL_BUILDDIR) 
	cp $(SRCDIR)/Makefile $(OPENCL_BUILDDIR)/Makefile
	
$(OPENCL_BUILDDIR)/libhuayno_cl.a: $(OPENCL_BUILDDIR)/Makefile $(LIBFILES) src/Makefile
	make -C $(OPENCL_BUILDDIR) libhuayno_cl.a VPATH=../src CFLAGS="$(CFLAGS)" CXXFLAGS="$(CXXFLAGS)" LIBS="$(LIBS)" INCLUDE="$(INCLUDE)"
 
	
worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py HuaynoInterface -o $@
	
worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py HuaynoInterface -o $@

huayno_worker: __init__.py worker_code.cc worker_code.h  $(BUILDDIR)/libhuayno.a $(BUILDDIR)/interface.o
	$(MPICXX) $(SC_FLAGS) $(CXXFLAGS) $(INCLUDE) worker_code.cc $(BUILDDIR)/interface.o $(BUILDDIR)/libhuayno.a -o $@ $(LIBS) $(SC_CLIBS)
	
huayno_worker_mp: CFLAGS   += $(OPENMP_CFLAGS) 
huayno_worker_mp: CXXFLAGS   += $(OPENMP_CFLAGS)	
huayno_worker_mp: __init__.py worker_code.cc worker_code.h $(OPENMP_BUILDDIR)/libhuayno.a $(OPENMP_BUILDDIR)/interface.o 
	$(MPICXX) $(SC_FLAGS) $(CXXFLAGS) $(INCLUDE) worker_code.cc $(OPENMP_BUILDDIR)/interface.o $(OPENMP_BUILDDIR)/libhuayno.a -o $@ $(LIBS) $(SC_CLIBS)
	

huayno_worker_cl: CFLAGS   += -fopenmp -DEVOLVE_OPENCL
huayno_worker_cl: CXXFLAGS   += -fopenmp -DEVOLVE_OPENCL
huayno_worker_cl: LIBS += $(CL_LIBS)  -lstdcl
huayno_worker_cl: INCLUDE += $(CL_FLAGS)
huayno_worker_cl: __init__.py worker_code.cc worker_code.h $(OPENCL_BUILDDIR)/libhuayno_cl.a $(OPENCL_BUILDDIR)/interface.o 
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) worker_code.cc $(OPENCL_BUILDDIR)/interface.o $(OPENCL_BUILDDIR)/libhuayno_cl.a -o $@ $(LIBS) $(SC_CLIBS)

.cc.o: $<
	$(MPICXX) -I$(SRCDIR)  $(SC_FLAGS) $(CXXFLAGS) $(INCLUDE) -c -o $@ $< 

$(BUILDDIR)/%.o: %.c
	$(MPICC) -I$(SRCDIR)  $(SC_FLAGS) $(CFLAGS) $(INCLUDE) -c -o $@ $< 
	
$(OPENMP_BUILDDIR)/%.o: %.c
	$(MPICC) -I$(SRCDIR)  $(SC_FLAGS) $(CFLAGS) $(INCLUDE) -c -o $@ $< 
	

$(OPENCL_BUILDDIR)/%.o: %.c
	$(MPICC) -I$(SRCDIR)  $(SC_FLAGS) $(CFLAGS) $(INCLUDE) -c -o $@ $< 
	
.c.o: $<
	$(MPICC) -I$(SRCDIR)  $(SC_FLAGS) $(CFLAGS) $(INCLUDE) -c -o $@ $< 

%.clh: %.cl
	awk 'BEGIN{print "const char srcstr[]=" } {if(substr($$0,length($0))=="\\"){$$0=$$0"\\"};print "\""$$0"\\n\""} END{print ";"}' $< > $@
