# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=./././.
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g -DCOLLISIONS_DUMMY -DLIBREBOUND 
CXXFLAGS += $(CFLAGS)  
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/librebound.a
CODEDIR = src/rebound
LIB_OBJ = $(CODEDIR)/src/collision.o $(CODEDIR)/src/particle.o $(CODEDIR)/src/tools.o $(CODEDIR)/src/integrator.o $(CODEDIR)/src/integrator_ias15.o $(CODEDIR)/src/integrator_whfast.o $(CODEDIR)/src/integrator_leapfrog.o $(CODEDIR)/src/integrator_sei.o $(CODEDIR)/src/integrator_hybrid.o $(CODEDIR)/src/integrator_wh.o $(CODEDIR)/src/gravity.o $(CODEDIR)/src/rebound.o $(CODEDIR)/src/input.o $(CODEDIR)/src/output.o $(CODEDIR)/src/boundary.o $(CODEDIR)/src/tree.o

DOWNLOAD_FROM_WEB = $(PYTHON) ./download_http.py
PATCH_FILES = $(PYTHON) ./patch_files.py
CODE_GENERATOR = $(AMUSE_DIR)/build.py

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
	@echo "DOWNLOAD_CODES is not set. Rebound will not be downloaded and build."
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
	$(RM) -f rebound_cython.* rebound_cython_worker rebound_cython_worker.*
	$(RM) -f rebound_worker_cython
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
	$(RM) rebound_cython.*

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
	

rebound_worker_cython: rebound_cython.so
	$(CODE_GENERATOR) --type=cython -m script -x amuse.community.rebound.interface ReboundInterface -o $@ --cython-import rebound_cython
	
rebound_cython.so: rebound_cython.o $(CODELIB) $(OBJS)
	$(MPICXX) -shared $(CXXFLAGS) $(PYTHONDEV_LDFLAGS) $(AM_CFLAGS) $(SC_FLAGS) $(LDFLAGS) rebound_cython.o  -o $@ $(SC_CLIBS) $(AM_LIBS) $(OBJS) $(CODELIB)

rebound_cython.o: rebound_cython.cc
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(AM_CFLAGS) $(PYTHONDEV_CFLAGS) -c -o $@ $< 

rebound_cython.pyx: interface.py worker_code.h
	$(CODE_GENERATOR) --type=cython -m mpi interface.py ReboundInterface -o $@

rebound_cython.cc: rebound_cython.pyx
	$(CYTHON) $< -o $@



.cc.o: $<
	$(CXX) $(CXXFLAGS) $(SC_FLAGS) -I$(CODEDIR)/src -c -o $@ $< 
