MPICXX ?= mpicxx
CXXFLAGS += -Wall -g -DTOOLBOX

AMUSE_DIR?=../../../..
CODE_GENERATOR = $(AMUSE_DIR)/build.py

SRCMAKEDIR   = src/mmas2
SRCDIR       = src/mmas2/src
MMAS_SRCDIR  = src/mmas2/src/mmas
EOS_SRCDIR   = src/mmas2/src/eos
USM_SRCDIR   = src/mmas2/src/usm
STD_SRCDIR   = src/mmas2/src/std

#INCLUDES = -I/opt/local/include -Isrc/mmas2/src
INCLUDES = -Isrc/mmas2/src
LIBS = -lgsl -lgslcblas \
	-Lsrc/mmas2/src/mmas -l_mmas2 \
	-Lsrc/mmas2/src/eos -leos \
	-Lsrc/mmas2/src/usm -lusm \
	-Lsrc/mmas2/src/std -lstd

OBJECTS = interface.o
SOURCES = interface.cc


all: mmas2 mmams_worker

mmas2:
	cd $(SRCMAKEDIR) && cmake . && make 

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
	make -C $(SRCMAKEDIR) clean
	rm -f *.pyc *.o
	rm -f worker_code.cc worker_code.h mmams_worker
	rm -f worker_code-sockets.cc mmams_worker_sockets

distclean:       clean
	rm -fr ccache
	cd $(SRCMAKEDIR)  && rm -rf CMakeCache.txt Makefile *.cmake *~ CMakeFiles
	cd $(SRCDIR)      && rm -rf CMakeCache.txt Makefile *.cmake *~ CMakeFiles
	cd $(MMAS_SRCDIR) && rm -rf CMakeCache.txt Makefile *.cmake *~ CMakeFiles ccache
	cd $(EOS_SRCDIR)  && rm -rf CMakeCache.txt Makefile *.cmake *~ CMakeFiles ccache
	cd $(USM_SRCDIR)  && rm -rf CMakeCache.txt Makefile *.cmake *~ CMakeFiles ccache
	cd $(STD_SRCDIR)  && rm -rf CMakeCache.txt Makefile *.cmake *~ CMakeFiles ccache
