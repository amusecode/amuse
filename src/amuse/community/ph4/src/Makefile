
# Makefile to build the libph4.a library.  The first draft was
# auto-generated; subsequent versions are hand-coded.  Save a copy!

# Optimization may cause failure on some systems.  Replace by -O0 in
# that case.

OPT	  = -O4

GPU	  = -DGPU
CXX	  = openmpicxx
CFLAGS   += -g -Wall $(OPT)
CXXFLAGS += $(CFLAGS) 

CODELIB = libph4.a

CODELIB_GPU = libph4gpu.a

CODEOBJS = debug.o diag.o gpu.o grape.o idata.o jdata.o scheduler.o \
	close_encounter.o two_body.o util.o

CODEOBJS_GPU = $(subst .o,.gpuo, $(CODEOBJS))

GPUFLAGS = $(GPU)

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB) 

gpulib: $(CODELIB_GPU)

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

$(CODELIB_GPU): $(CODEOBJS_GPU)
	echo $(CODEOBJS_GPU)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS_GPU)
	$(RANLIB) $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

$(CODEOBJS_GPU): %.gpuo: %.cc
	$(CXX) $(CXXFLAGS) $(GPUFLAGS) -c -o $@ $< 

tar:
	tar czf ph4.tgz Makefile.ph4 *.h *.cc plummer*.in*

clean:
	$(RM) -r -f *.o *~ *.tgz *.a *.gpuo
