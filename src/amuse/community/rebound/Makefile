# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g -DCOLLISIONS_DUMMY -DLIBREBOUND
CXXFLAGS += $(CFLAGS)  
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

ifdef OPENMP_ENABLED
export OPENMP=1
CFLAGS   += $(OPENMP_CFLAGS)
endif

OBJS = interface.o

CODELIB = src/librebound.a
CODEDIR = src/rebound
LIB_OBJ = \
	  $(CODEDIR)/src/rebound.o \
	  $(CODEDIR)/src/tree.o \
	  $(CODEDIR)/src/particle.o \
	  $(CODEDIR)/src/gravity.o \
          $(CODEDIR)/src/integrator.o \
	  $(CODEDIR)/src/integrator_whfast.o \
	  $(CODEDIR)/src/integrator_whfast512.o \
	  $(CODEDIR)/src/integrator_saba.o \
          $(CODEDIR)/src/integrator_ias15.o \
	  $(CODEDIR)/src/integrator_sei.o \
	  $(CODEDIR)/src/integrator_bs.o \
	  $(CODEDIR)/src/integrator_leapfrog.o \
          $(CODEDIR)/src/integrator_mercurius.o \
          $(CODEDIR)/src/integrator_eos.o \
	  $(CODEDIR)/src/boundary.o \
	  $(CODEDIR)/src/input.o \
	  $(CODEDIR)/src/binarydiff.o \
	  $(CODEDIR)/src/output.o \
          $(CODEDIR)/src/collision.o \
          $(CODEDIR)/src/communication_mpi.o \
	  $(CODEDIR)/src/display.o \
	  $(CODEDIR)/src/tools.o \
	  $(CODEDIR)/src/rotations.o \
	  $(CODEDIR)/src/derivatives.o \
          $(CODEDIR)/src/simulationarchive.o \
          $(CODEDIR)/src/glad.o \
	  $(CODEDIR)/src/integrator_janus.o \
	  $(CODEDIR)/src/transformations.o \
	  $(CODEDIR)/src/fmemopen.o \
          $(CODEDIR)/src/server.o

DOWNLOAD_FROM_WEB = $(PYTHON) ./download.py
PATCH_FILES = $(PYTHON) ./patch_files.py

AM_LIBS = -L$(AMUSE_DIR)/lib/amuse_mpi -lamuse_mpi
AM_CFLAGS = -I$(AMUSE_DIR)/lib/amuse_mpi


all: rebound_worker 

$(CODEDIR)/Makefile:
	make -C . download

download:
	$(RM) -Rf .pc
	$(RM) -Rf src
	mkdir src
	$(DOWNLOAD_FROM_WEB)
	$(PATCH_FILES)

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ rebound_worker worker_code.cc
	$(RM) -Rf $(CODEDIR)/src/*.o
	$(RM) -Rf $(CODEDIR)/src/*.so
	$(RM) -Rf src/librebound.a
	-make -C $(CODEDIR) clean
	
distclean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h *.mod
	$(RM) *~ rebound_worker
	$(RM) -Rf src
	$(RM) -Rf src.*
	$(RM) -Rf $(CODEDIR)/src/*.o
	$(RM) -Rf $(CODEDIR)/src/*.so
	$(RM) -Rf src/librebound.a
	$(RM) -Rf .pc

$(CODELIB): $(CODEDIR)/Makefile
	make -C $(CODEDIR) librebound CC=$(CC)
	ar rv $(CODELIB) $(LIB_OBJ)
	ranlib $(CODELIB)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py ReboundInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py  -i amuse.community.interface.stopping_conditions.StoppingConditionInterface ReboundInterface -o $@

rebound_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(LDFLAGS) -I$(CODEDIR)/src $< $(OBJS) $(CODELIB) $(SC_CLIBS)  -o $@
	
.cc.o: $<
	$(CXX) $(CXXFLAGS) $(SC_FLAGS) -I$(CODEDIR)/src -c -o $@ $< 
