# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx

CODE_GENERATOR = $(AMUSE_DIR)/build.py

BUILDDIR = build
SRCMAKEDIR   = src/mmas2
SRCDIR       = src/mmas2/src
MMAS_DIR  = src/mmas2/src/mmas
EOS_DIR   = src/mmas2/src/eos
USM_DIR   = src/mmas2/src/usm
STD_SRCMMAS_DIRDIR   = src/mmas2/src/std
GSL_LIBS ?= -lgsl -lgslcblas -lm

#INCLUDES = -I/opt/local/include -Isrc/mmas2/src
INCLUDES = $(GSL_FLAGS)  -Isrc/mmas2/src
LIBS =  -L$(BUILDDIR)/src/mmas -l_mmas2 \
	-L$(BUILDDIR)/src/eos -leos \
	-L$(BUILDDIR)/src/usm -lusm \
	-L$(BUILDDIR)/src/std -lstd \
        $(GSL_LIBS)

OBJECTS = interface.o
SOURCES = interface.cc

CMAKE_GENERATOR ?= Unix Makefiles

all: $(BUILDDIR)/Makefile  mmas2 mmams_worker

$(BUILDDIR)/Makefile:
	@-mkdir build
	cd $(BUILDDIR) && cmake ../src/mmas2 -G "$(CMAKE_GENERATOR)" -DGSL_FLAGS=$(GSL_FLAGS)
	
mmas2: $(BUILDDIR)/Makefile
	make -C $(BUILDDIR)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c $^ MakeMeAMassiveStarInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H $^ MakeMeAMassiveStarInterface -o $@

mmams_worker: worker_code.cc worker_code.h $(OBJECTS)
	$(MPICXX) $(CXXFLAGS) $< $(OBJECTS) -o $@  $(LIBS)

interface.o: interface.cc
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $< $(LIBS)

clean:
	@-make -C $(SRCMAKEDIR) clean
	cd $(SRCMAKEDIR)  && rm -rf CMakeCache.txt Makefile *.cmake *~ CMakeFiles
	cd $(SRCDIR)      && rm -rf CMakeCache.txt Makefile *.cmake *~ CMakeFiles
	cd $(MMAS_DIR) && rm -rf CMakeCache.txt Makefile *.cmake *~ CMakeFiles ccache
	cd $(EOS_DIR)  && rm -rf CMakeCache.txt Makefile *.cmake *~ CMakeFiles ccache
	cd $(USM_DIR)  && rm -rf CMakeCache.txt Makefile *.cmake *~ CMakeFiles ccache
	cd $(STD_DIR)  && rm -rf CMakeCache.txt Makefile *.cmake *~ CMakeFiles ccache
	rm -f *.pyc *.o
	rm -f worker_code.cc worker_code.h mmams_worker
	rm -f worker_code-sockets.cc mmams_worker_sockets
	rm -rf build

distclean: clean

