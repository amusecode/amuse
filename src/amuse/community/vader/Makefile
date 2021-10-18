# standard amuse configuration include
# config.mk will be made after ./configure has run
#AMUSE_DIR?=../../../..
AMUSE_DIR?=$(shell amusifier --get-amuse-dir)
#AMUSE_DIR?= $(AMUSE_DIR)
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx
MPICC    ?= mpicc

CFLAGS   += -Wall -g 
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm -lgsl -lgslcblas $(MUSE_LD_FLAGS) 

TARGETS := $(patsubst src/prob/userFunc_%.c, vader_worker_%, ${shell ls src/prob/userFunc_*.c})

OBJS = interface.o

CODELIB = src/libvader_%.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

PROB = none

BUILDDIR = build

#all: vader_worker_% 
all: ${TARGETS}

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ vader_worker_* worker_code.cc
	$(RM) -Rf build_*
	make -C src clean

vader_worker: ${TARGETS}

${TARGETS}: vader_worker_%: worker_code.cc worker_code.h interface_%.o $(BUILDDIR)_%/libvader.a
	$(MPICXX) $(CXXFLAGS) $< interface_$*.o $(BUILDDIR)_$*/libvader.a -o $@ $(LDFLAGS)

$(BUILDDIR)_%:
	-mkdir $@

$(BUILDDIR)_%/Makefile: src/Makefile | $(BUILDDIR)_%
	cp src/Makefile $@

$(BUILDDIR)_%/userFunc.h: src/userFunc.h | $(BUILDDIR)_%
	cp src/userFunc.h $@

$(BUILDDIR)_%/vader_common.h: src/vader_common.h | $(BUILDDIR)_%
	cp src/vader_common.h $@

$(BUILDDIR)_%/libvader.a: $(BUILDDIR)_%/Makefile $(BUILDDIR)_%/userFunc.h $(BUILDDIR)_%/vader_common.h
	make -C $(BUILDDIR)_$* all PROB=$* VPATH=../src

#$(CODELIB):
#	make -C src all

interface_%.o: interface.cc worker_code.h
	$(MPICC) $(CXXFLAGS) -c -o $@ $<

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py VaderInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py VaderInterface -o $@

#vader_worker_%: worker_code_%.cc worker_code_%.h $(CODELIB) $(OBJS)
#	$(MPICXX) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@ $(LDFLAGS) 

.cc.o: $<
	$(MPICC) $(CXXFLAGS) -c -o $@ $< 
