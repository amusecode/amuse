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
LIB_OBJ = $(CODEDIR)/src/rebound.o $(CODEDIR)/src/tree.o $(CODEDIR)/src/particle.o $(CODEDIR)/src/gravity.o \
          $(CODEDIR)/src/integrator.o $(CODEDIR)/src/integrator_whfast.o $(CODEDIR)/src/integrator_whfasthelio.o \
          $(CODEDIR)/src/integrator_ias15.o $(CODEDIR)/src/integrator_sei.o $(CODEDIR)/src/integrator_leapfrog.o \
          $(CODEDIR)/src/integrator_hermes.o $(CODEDIR)/src/boundary.o \
          $(CODEDIR)/src/collision.o $(CODEDIR)/src/tools.o \
          $(CODEDIR)/src/communication_mpi.o $(CODEDIR)/src/display.o $(CODEDIR)/src/derivatives.o \
          $(CODEDIR)/src/glad.o $(CODEDIR)/src/integrator_janus.o $(CODEDIR)/src/transformations.o \
          $(CODEDIR)/src/simulationarchive.o $(CODEDIR)/src/output.o \
          $(CODEDIR)/src/input.o \

DOWNLOAD_FROM_WEB = $(PYTHON) ./download.py
PATCH_FILES = $(PYTHON) ./patch_files.py
CODE_GENERATOR ?= $(AMUSE_DIR)/build.py

AM_LIBS = -L$(AMUSE_DIR)/lib/amuse_mpi -lamuse_mpi
AM_CFLAGS = -I$(AMUSE_DIR)/lib/amuse_mpi


all: rebound_worker 


ifdef DOWNLOAD_CODES
$(CODEDIR)/Makefile:
	make -C . download
else
$(CODEDIR)/Makefile:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. Rebound will not be downloaded and built."
	@echo "If you do want Rebound, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif

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
	make -C $(CODEDIR) all CC=$(CC)
	ar rv $(CODELIB) $(LIB_OBJ)
	ranlib $(CODELIB)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py ReboundInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py  -i amuse.support.codes.stopping_conditions.StoppingConditionInterface ReboundInterface -o $@

rebound_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(LDFLAGS) -I$(CODEDIR)/src $< $(OBJS) $(CODELIB) $(SC_CLIBS)  -o $@
	
.cc.o: $<
	$(CXX) $(CXXFLAGS) $(SC_FLAGS) -I$(CODEDIR)/src -c -o $@ $< 
