
# Makefile to build the libsmallN.a library.  The first draft was
# auto-generated; subsequent versions are hand-coded.  Save a copy!

# Optimization may cause failure on some systems.  Replace by -O0 in
# that case.

AMUSE_DIR ?= ../../../../..
OPT	 = -O4 -I$(AMUSE_DIR)/lib/stopcond

CXX	 ?= mpicxx
CFLAGS   += -g -Wall $(OPT)
CXXFLAGS += $(CFLAGS) 

CODELIB = libsmallN.a

CODEOBJS = smallN-lib.o smallN_unpert.o kepler-lib.o hdyn.o \
	   analyze-lib.o nstab.co f2c.co util.o \
	   smallN_interface.o

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

#-----------------------------------------------------------------------
#
# Special cases:

nstab.co:    nstab.c
	$(MPICC) $(CFLAGS) -c nstab.c -o nstab.co

f2c.co:    f2c.c
	$(MPICC) $(CFLAGS) -c f2c.c -o f2c.co

smallN-lib.o:	smallN.cc
	$(MPICXX) $(CFLAGS) -c smallN.cc -o smallN-lib.o

kepler-lib.o:	kepler.cc
	$(MPICXX) $(CFLAGS) -c kepler.cc -o kepler-lib.o

analyze-lib.o:	analyze.cc
	$(MPICXX) $(CFLAGS) -c analyze.cc -o analyze-lib.o

#-----------------------------------------------------------------------

tar:
	tar czf smallN.tgz Makefile *.h *.cc plummer*.in*

clean:
	$(RM) -r -f *.o *.co *~ *.tgz *.a
