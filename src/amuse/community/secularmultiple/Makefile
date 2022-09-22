ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif

-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx
MPICC ?= mpicc

CC = $(MPICC)
CXX = $(MPICXX)

CXXFLAGS += -g -DTOOLBOX -O3 $(MUSE_INCLUDE_DIR)
LDFLAGS+ = -lm $(MUSE_LD_FLAGS)

SRCDIR=src

OBJS = $(SRCDIR)/interface.o $(SRCDIR)/src/types.o $(SRCDIR)/src/evolve.o  \
$(SRCDIR)/src/structure.o $(SRCDIR)/src/ODE_system.o \
$(SRCDIR)/src/root_finding.o $(SRCDIR)/src/newtonian.o \
$(SRCDIR)/src/postnewtonian.o $(SRCDIR)/src/tides.o \
$(SRCDIR)/src/external.o $(SRCDIR)/src/VRR.o $(SRCDIR)/src/cvode/cvode.o \
$(SRCDIR)/src/cvode/cvode_dense.o $(SRCDIR)/src/cvode/cvode_direct.o \
$(SRCDIR)/src/cvode/cvode_io.o $(SRCDIR)/src/cvode/nvector_serial.o \
$(SRCDIR)/src/cvode/sundials_dense.o $(SRCDIR)/src/cvode/sundials_direct.o \
$(SRCDIR)/src/cvode/sundials_math.o $(SRCDIR)/src/cvode/sundials_nvector.o 

all: secularmultiple_worker

clean:
	rm -f *.so *.o *.pyc worker_code.cc $(SRCDIR)/*.o* $(SRCDIR)/src/*.o*
	
distclean: clean
	rm -rf src
	
worker_code.cc: worker_code.h interface.py
	$(CODE_GENERATOR) --type=c interface.py SecularMultipleInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py SecularMultipleInterface -o $@

secularmultiple_worker: worker_code.cc $(OBJS)
	$(MPICXX) $(CXXFLAGS) worker_code.cc $(OBJS) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.c.o: $<
	$(CC) $(CXXFLAGS) -c -o $@ $<

src/interface.cpp:
	make -C . download

download:
	git clone -b master https://github.com/hamers/secularmultiple.git src
