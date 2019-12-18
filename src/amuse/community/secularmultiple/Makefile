ifndef AMUSE_DIR
    AMUSE_DIR=$(shell amusifier --get-amuse-dir)
endif

-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx
MPICC ?= mpicc

CC = $(MPICC)
CXX = $(MPICXX)

SRCDIR=secularmultiple

CODE_GENERATOR = $(AMUSE_DIR)/build.py

CXXFLAGS += -g -DTOOLBOX -O3 $(MUSE_INCLUDE_DIR)
LDFLAGS+ = -lm $(MUSE_LD_FLAGS)

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

cleanall: clean
	$(RM) worker_code *

clean:
	rm -f *.so *.o *.pyc worker_code.cc $(SRCDIR)/src/*.o*
	
distclean: clean
	rm -rf $(SRCDIR)
	
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

ifdef DOWNLOAD_CODES
secularmultiple/interface.cpp:
	make -C . download
else
secularmultiple/interface.cpp:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. SecularMultiple will not be downloaded and build."
	@echo "If you do want it, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif 

download:
	git clone -b master https://github.com/hamers/secularmultiple.git $(SRCDIR)
