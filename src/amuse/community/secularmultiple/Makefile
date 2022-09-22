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

OBJS = interface.o ./src/types.o ./src/evolve.o  \
./src/structure.o ./src/ODE_system.o \
./src/root_finding.o ./src/newtonian.o \
./src/postnewtonian.o ./src/tides.o \
./src/external.o ./src/VRR.o ./src/cvode/cvode.o \
./src/cvode/cvode_dense.o ./src/cvode/cvode_direct.o \
./src/cvode/cvode_io.o ./src/cvode/nvector_serial.o \
./src/cvode/sundials_dense.o ./src/cvode/sundials_direct.o \
./src/cvode/sundials_math.o ./src/cvode/sundials_nvector.o 

all: secularmultiple_worker

clean:
	rm -f *.so *.o *.pyc worker_code.cc ./src/*.o*
	
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
	git clone -b master https://github.com/hamers/secularmultiple.git src
