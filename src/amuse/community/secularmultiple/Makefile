ifndef AMUSE_DIR
    AMUSE_DIR=../../../..
endif

-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx
MPICC ?= mpicc

CC = $(MPICC)
CXX = $(MPICXX)

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CXXFLAGS += -g -DTOOLBOX -O3 $(MUSE_INCLUDE_DIR)
LDFLAGS+ = -lm $(MUSE_LD_FLAGS)

OBJS = interface.o src/types.o src/evolve.o src/structure.o src/ODE_system.o src/root_finding.o src/newtonian.o src/postnewtonian.o src/tides.o src/external.o src/cvode/cvode.o src/cvode/cvode_dense.o src/cvode/cvode_direct.o src/cvode/cvode_io.o src/cvode/nvector_serial.o src/cvode/sundials_dense.o src/cvode/sundials_direct.o src/cvode/sundials_math.o src/cvode/sundials_nvector.o 

all: worker_code

cleanall: clean
	$(RM) worker_code *

clean:
	rm -f *.so *.o *.pyc worker_code.cc src/*.o*
	
distclean:
	rm -f *.so *.o *.pyc worker_code.cc src/*.o* src/cvode/*.o*

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py SecularMultipleInterface -o $@

worker_code: worker_code.cc $(OBJS)
	$(MPICXX) $(CXXFLAGS) $@.cc $(OBJS) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.c.o: $<
	$(CC) $(CXXFLAGS) -c -o $@ $<
