
# Makefile to build the libph4.a library.  The first draft was
# auto-generated; subsequent versions are hand-coded.  Save a copy!

# Optimizaton may cause failure on some systems.  Replace by -O0 in
# that case.

OPT	  = -O4

GPU	  = -DGPU
CXX	  = openmpicxx
CFLAGS   += -Wall -g $(GPU) $(OPT)
CXXFLAGS += $(CFLAGS) 

CODELIB = libph4.a

CODEOBJS = diag.o gpu.o grape.o idata.o jdata.o scheduler.o

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB) 

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

tar:
	tar czf ph4.tgz Makefile.ph4 *.h *.cc plummer*.in*

clean:
	$(RM) -r -f *.o *~ *.tgz *.a
