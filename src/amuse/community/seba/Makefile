# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
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

CODE_GENERATOR ?= $(AMUSE_DIR)/build.py

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

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
