
# Makefile to build the libph4.a library.  The first draft was
# auto-generated; subsequent versions are hand-coded.  Save a copy!

# Optimization may cause failure on some systems.  Replace by -O0 in
# that case.

AMUSE_DIR?=../../../../..
OPT	  = -O4 -I$(AMUSE_DIR)/lib/stopcond

GPU	  = -DGPU
MPICXX    ?= mpicxx
MPICC     ?= mpicc
CXX       ?= mpicxx
CFLAGS   += -g -Wall $(OPT)
CXXFLAGS += $(CFLAGS) 

CODELIB = libph4.a

CODELIB_GPU = libph4gpu.a

CODEOBJS = debug.o diag.o gpu.o grape.o idata.o jdata.o scheduler.o \
	   close_encounter.o two_body.o util.o
CODEOBJS_GPU = $(subst .o,.gpuo, $(CODEOBJS))

# Needed if we don't use AMUSE stopping conditions:

EXTRAOBJS = kepler-lib.o hdyn.o smallN-lib.o smallN_unpert.o \
	    analyze-lib.o nstab.co f2c.co

GPUFLAGS = $(GPU)

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB) 

gpulib: $(CODELIB_GPU)

$(CODELIB): $(CODEOBJS) $(EXTRAOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS) $(EXTRAOBJS)
	$(RANLIB) $@

$(CODELIB_GPU): $(CODEOBJS_GPU) $(EXTRAOBJS)
	echo $(CODEOBJS_GPU)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS_GPU) $(EXTRAOBJS)
	$(RANLIB) $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

$(CODEOBJS_GPU): %.gpuo: %.cc
	$(CXX) $(CXXFLAGS) $(GPUFLAGS) -c -o $@ $< 

#-----------------------------------------------------------------------

# Should only need to build nstab.c after a (rare) update of nstab.f
# from R. Mardling.
update_nstab:
	f2c nstab.f

nstab.co:    nstab.c
	$(MPICC) $(CFLAGS) -c nstab.c -o nstab.co

f2c.co:    f2c.c
	$(MPICC) $(CFLAGS) -c f2c.c -o f2c.co

#-----------------------------------------------------------------------

smallN-lib.o:	smallN.cc
	$(MPICXX) $(CFLAGS) -c smallN.cc -o smallN-lib.o

kepler-lib.o:	kepler.cc
	$(MPICXX) $(CFLAGS) -c kepler.cc -o kepler-lib.o

analyze-lib.o:	analyze.cc
	$(MPICXX) $(CFLAGS) -c analyze.cc -o analyze-lib.o

#-----------------------------------------------------------------------

tar:
	tar czf ph4.tgz Makefile.ph4 *.h *.cc plummer*.in*

clean:
	$(RM) -r -f *.o *.co *~ *.tgz *.a *.gpuo
