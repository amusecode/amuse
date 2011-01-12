MPICXX ?= mpicxx
MPICC ?= mpicc

VERSION = 4.0

ATHENA_DIR = src/athena$(VERSION)

CC = $(MPICC)
CFLAGS   += -g -I$(ATHENA_DIR)/src 
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -ldl -lm

AMUSE_OBJS = interface_$(VERSION).o
AMUSE_OBJS_SELFGRAV = interface_$(VERSION).selfgrav_o

CODELIB = $(ATHENA_DIR)/bin/libathena.a

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py
DOWNLOAD_FROM_WEB = ./downloadhttp.py
PATCH_FILES = ./patch_files.py

SCLIBS ?= -L../../../../lib/stopcond -lstopcond
SCINC ?= -I../../../../lib/stopcond

ATHENA_OBJS = $(subst $(ATHENA_DIR)/src/main.o,,$(wildcard $(ATHENA_DIR)/src/*.o $(ATHENA_DIR)/src/*/*.o))

FFTWLIB = -L/usr/lib -lfftw3
FFTWINC = -I/usr/include
  
all: athena_worker 

$(ATHENA_DIR):
ifdef DOWNLOAD_CODES
	$(RM) -Rf .pc
	$(DOWNLOAD_FROM_WEB) $(VERSION)
	#cd $(ATHENA_DIR)/src;patch Makefile.in < ../../../makefile.patch
	$(PATCH_FILES)
else
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


clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) -f *~ athena_worker
	$(RM) -f $(CODELIB) $(AMUSE_OBJS_SELFGRAV)

codeclean: clean
	make -C $(ATHENA_DIR) clean

distclean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code
	$(RM) -Rf src

$(CODELIB): $(ATHENA_DIR) 
	cd $(ATHENA_DIR); ./configure --with-gas=hydro --enable-mpi --enable-smr
	cd $(ATHENA_DIR); ln -fs ../../../amuse_problem_$(VERSION).c src/problem.c
	cd $(ATHENA_DIR); make all CFLAGS=-g CC=$(MPICC) LDR=$(MPICC)
	rm -f $(CODELIB)
	cd $(ATHENA_DIR); make all CFLAGS=-g CC=$(MPICC) LDR=$(MPICC)
	

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py AthenaInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py -i amuse.support.codes.stopping_conditions.StoppingConditionInterface AthenaInterface -o $@

athena_worker: worker_code.cc worker_code.h $(CODELIB) $(AMUSE_OBJS)
	$(MPICXX) $(CXXFLAGS) $(SCINC) worker_code.cc $(AMUSE_OBJS) $(CODELIB) $(LDFLAGS) -o $@ $(SCLIBS)

athena_worker_selfgrav: worker_code.cc worker_code.h interface_$(VERSION).c
	cd $(ATHENA_DIR); make clean
	$(RM) -f $(CODELIB)
	cd $(ATHENA_DIR); ./configure --with-gas=hydro --with-gravity=fft --enable-fft

	
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
	
.c.o: $<
	$(CC) $(CXXFLAGS) $(SCINC) -c -o $@ $< 
	


