# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx

CFLAGS   += -Isrc/include -Isrc/include/star $(SC_FLAGS)
CXXFLAGS += $(CFLAGS)
LDFLAGS  += -Lsrc/node  -Lsrc/std 
LIBS = -lnode -lstd -lm
# so that seba build process works
unexport CFLAGS CXXFLAGS LDFLAGS

AMUSE_OBJS = interface.o

AM_LIBS = -L$(AMUSE_DIR)/lib/amuse_mpi -lamuse_mpi
AM_CFLAGS = -I$(AMUSE_DIR)/lib/amuse_mpi


CODELIB =  src/sstar/libsstar.a src/dstar/libdstar.a 

CODE_GENERATOR = $(AMUSE_DIR)/build.py
DOWNLOAD_FROM_SVN = ./downloadsvn.py

all: seba_worker

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) *~ worker_code seba_worker
	@-make -C src clean

distclean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) *~ worker_code

$(CODELIB): src
	cd src; make -C. CXX="$(CXX) $(CXXFLAGS)"

worker_code.cc: interface.py
	CFLAGS="$(CFLAGS)" $(CODE_GENERATOR) --type=c interface.py SeBaInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py SeBaInterface -o $@

seba_worker: worker_code.cc worker_code.h $(CODELIB) $(AMUSE_OBJS)
	$(MPICXX) $(CXXFLAGS)  $(LDFLAGS) $< $(AMUSE_OBJS) $(CODELIB) $(LIBS) -o $@  $(LIBS) $(SC_CLIBS) 


seba_worker_cython: seba_cython.so
	$(CODE_GENERATOR) --type=cython -m script -x amuse.community.seba.interface SeBaInterface -o $@ --cython-import seba_cython
	
seba_cython.so: seba_cython.o $(CODELIB) $(AMUSE_OBJS)
	$(MPICXX) -shared $(CXXFLAGS) $(PYTHONDEV_LDFLAGS) $(AM_CFLAGS) $(SC_FLAGS) $(LDFLAGS) $(CODELIB) seba_cython.o  -o $@ $(SC_CLIBS) $(AM_LIBS) $(AMUSE_OBJS) $(CODELIB) $(LIBS)

seba_cython.o: seba_cython.cc
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(AM_CFLAGS) $(PYTHONDEV_CFLAGS) -c -o $@ $< 

seba_cython.pyx: interface.py worker_code.h
	$(CODE_GENERATOR) --type=cython -m mpi interface.py SeBaInterface -o $@

seba_cython.cc: seba_cython.pyx
	$(CYTHON) --cplus $< -o $@


.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
