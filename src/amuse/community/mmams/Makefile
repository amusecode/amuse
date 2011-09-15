MPICXX ?= mpicxx

AMUSE_DIR?=../../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py

BUILDDIR = build
SRCMAKEDIR   = src/mmas2
SRCDIR       = src/mmas2/src
MMAS_DIR  = src/mmas2/src/mmas
EOS_DIR   = src/mmas2/src/eos
USM_DIR   = src/mmas2/src/usm
STD_SRCMMAS_DIRDIR   = src/mmas2/src/std

#INCLUDES = -I/opt/local/include -Isrc/mmas2/src
INCLUDES = -Isrc/mmas2/src
LIBS = $(GSLFLAGS) -lgsl -lgslcblas \
	-L$(BUILDDIR)/src/mmas -l_mmas2 \
	-L$(BUILDDIR)/src/eos -leos \
	-L$(BUILDDIR)/src/usm -lusm \
	-L$(BUILDDIR)/src/std -lstd

OBJECTS = interface.o
SOURCES = interface.cc

all: $(BUILDDIR)/Makefile  mmas2 mmams_worker

$(BUILDDIR)/Makefile:
	@-mkdir build
	cd $(BUILDDIR) && cmake ../src/mmas2
	
mmas2: $(BUILDDIR)/Makefile
	make -C $(BUILDDIR)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c $^ MakeMeAMassiveStarInterface -o $@

worker_code-sockets.cc: interface.py
	$(CODE_GENERATOR) --type=c --mode=sockets $^ MakeMeAMassiveStarInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H $^ MakeMeAMassiveStarInterface -o $@

mmams_worker: worker_code.cc worker_code.h $(OBJECTS)
	$(MPICXX) $(CXXFLAGS) $< $(OBJECTS) -o $@  $(LIBS)

mmams_worker_sockets: worker_code-sockets.cc worker_code.h $(OBJECTS)
	$(MPICXX) $(CXXFLAGS) $< $(OBJECTS) -o $@  $(LIBS)

interface.o: interface.cc
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $< $(LIBS)

clean:
	@-make -C $(SRCMAKEDIR) clean
	rm -f *.pyc *.o
	rm -f worker_code.cc worker_code.h mmams_worker
	rm -f worker_code-sockets.cc mmams_worker_sockets
	rm -rf build

distclean: clean
	rm -rf ccache
