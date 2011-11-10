# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx
MPICC ?= mpicc
PYTHON ?= python

VERSION = 4.1

ATHENA_DIR = src/athena$(VERSION)

CC = $(MPICC)
CFLAGS   += -g -I$(ATHENA_DIR)/src 
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -ldl -lm

AMUSE_OBJS = interface_$(VERSION).o
AMUSE_OBJS_SELFGRAV = interface_$(VERSION).selfgrav_o
AMUSE_OBJS_MHD = interface_$(VERSION).mhd_o

CODELIB = $(ATHENA_DIR)/bin/libathena.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py
DOWNLOAD_FROM_WEB = $(PYTHON) ./downloadhttp.py
PATCH_FILES = $(PYTHON) ./patch_files.py

SCLIBS ?= -L$(AMUSE_DIR)/lib/stopcond -lstopcond
SCINC ?= -I$(AMUSE_DIR)/lib/stopcond

ATHENA_OBJS = $(subst $(ATHENA_DIR)/src/main.o,,$(wildcard $(ATHENA_DIR)/src/*.o $(ATHENA_DIR)/src/*/*.o))

FFTWLIB = -L/usr/lib -lfftw3
FFTWINC = -I/usr/include
 
.PHONY : download

all: athena_worker 

ifdef DOWNLOAD_CODES
$(ATHENA_DIR)/Makefile:
	make -C . download
else
$(ATHENA_DIR)/Makefile:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. Athena will not be downloaded and build."
	@echo "If you do want Athena, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif 

download:
	$(RM) -Rf .pc
	$(DOWNLOAD_FROM_WEB) $(VERSION)
	$(PATCH_FILES)

clean:
	$(RM) -f *.so *.o *.selfgrav_o *.pyc worker_code.cc worker_code.h 
	$(RM) -f *~ athena_worker  athena_worker_selfgrav
	$(RM) -f worker_code-sockets.cc athena_worker_sockets
	$(RM) -f $(CODELIB) $(AMUSE_OBJS_SELFGRAV)

codeclean: clean
	make -C $(ATHENA_DIR) clean

distclean: clean
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) -Rf .pc
	$(RM) *~ athena_worker
	$(RM) -Rf src

$(CODELIB): $(ATHENA_DIR)/Makefile
	cd $(ATHENA_DIR); ./configure --with-gas=hydro --enable-mpi --enable-smr 
	cd $(ATHENA_DIR); ln -fs ../../../amuse_problem_$(VERSION).c src/problem.c
	cd $(ATHENA_DIR); make all CFLAGS=-g CC=$(MPICC) LDR=$(MPICC)
	rm -f $(CODELIB)
	cd $(ATHENA_DIR); make all CFLAGS=-g CC=$(MPICC) LDR=$(MPICC)
	

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py AthenaInterface -o $@

worker_code-sockets.cc: interface.py
	$(CODE_GENERATOR) --type=c --mode=sockets interface.py AthenaInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py -i amuse.support.codes.stopping_conditions.StoppingConditionInterface AthenaInterface -o $@

athena_worker: worker_code.cc worker_code.h $(CODELIB) $(AMUSE_OBJS)
	$(MPICXX) $(CXXFLAGS) $(SCINC) worker_code.cc $(AMUSE_OBJS) $(CODELIB) $(LDFLAGS) -o $@ $(SCLIBS)

athena_worker_sockets: worker_code-sockets.cc worker_code.h $(CODELIB) $(AMUSE_OBJS)
	$(MPICXX) $(CXXFLAGS) $(SCINC) worker_code-sockets.cc $(AMUSE_OBJS) $(CODELIB) $(LDFLAGS) -o $@ $(SCLIBS)

athena_worker_selfgrav: worker_code.cc worker_code.h interface_$(VERSION).c
	cd $(ATHENA_DIR); make clean
	$(RM) -f $(CODELIB)
	cd $(ATHENA_DIR); ./configure --with-gas=hydro --with-gravity=fft --enable-fft --disable-cooling
	#cd $(ATHENA_DIR); ./configure --with-gas=hydro --with-gravity=fft --enable-fft --with-flux=roe --with-eos=adiabatic --with-order=2 --enable-h-correction --disable-cooling
	#cd $(ATHENA_DIR); ./configure --with-gas=hydro --with-gravity=multigrid

	
	$(CC) $(CXXFLAGS) $(SCINC) -c -o $(AMUSE_OBJS_SELFGRAV) interface_$(VERSION).c
	
	cd $(ATHENA_DIR); ln -fs ../../../amuse_problem_$(VERSION).c src/problem.c
	cd $(ATHENA_DIR); make all CFLAGS=-g CC=$(MPICC) LDR=$(MPICC)
	$(RM) -f $(CODELIB)
	cd $(ATHENA_DIR); make all CFLAGS=-g CC=$(MPICC) LDR=$(MPICC)
	$(MPICXX) $(CXXFLAGS) $(SCINC) worker_code.cc $(CODELIB) $(AMUSE_OBJS_SELFGRAV) $(CODELIB) $(LDFLAGS) -o $@ $(SCLIBS) $(FFTWLIB)
	$(RM) -f $(CODELIB)
	cd $(ATHENA_DIR); ./configure --with-gas=hydro --enable-mpi --enable-smr
	cd $(ATHENA_DIR); make clean
	$(RM) -f $(CODELIB)
	

athena_worker_mhd: worker_code.cc worker_code.h interface_$(VERSION).c
	cd $(ATHENA_DIR); make clean
	$(RM) -f $(CODELIB)
	cd $(ATHENA_DIR); ./configure --with-gas=mhd --with-eos=adiabatic --enable-mpi --enable-smr --with-order=3 --with-flux=hlld
	
	$(CC) $(CXXFLAGS) $(SCINC) -c -o $(AMUSE_OBJS_MHD) interface_$(VERSION).c
	
	cd $(ATHENA_DIR); ln -fs ../../../amuse_problem_$(VERSION).c src/problem.c
	cd $(ATHENA_DIR); make all CFLAGS=-g CC=$(MPICC) LDR=$(MPICC)
	$(RM) -f $(CODELIB)
	
	cd $(ATHENA_DIR); make all CFLAGS=-g CC=$(MPICC) LDR=$(MPICC)
	$(MPICXX) $(CXXFLAGS) $(SCINC) worker_code.cc $(CODELIB) $(AMUSE_OBJS_MHD) $(CODELIB) $(LDFLAGS) -o $@ $(SCLIBS) $(FFTWLIB)
	$(RM) -f $(CODELIB)
	cd $(ATHENA_DIR); ./configure --with-gas=hydro --enable-mpi --enable-smr
	cd $(ATHENA_DIR); make clean
	$(RM) -f $(CODELIB)
	

	
.c.o: $<
	$(CC) $(CXXFLAGS) $(SCINC) -c -o $@ $< 
	


