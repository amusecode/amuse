VERSION = 4.0

ATHENA_DIR = src/athena$(VERSION)

CC = mpicc
CFLAGS   += -g -I$(ATHENA_DIR)/src 
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -ldl -lm

AMUSE_OBJS = interface_$(VERSION).o

CODELIB = $(ATHENA_DIR)/bin/libathena.a

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py
DOWNLOAD_FROM_WEB = ./downloadhttp.py

SCLIBS ?= -L../../../../lib/stopcond -lstopcond
SCINC ?= -I../../../../lib/stopcond

ATHENA_OBJS = $(subst $(ATHENA_DIR)/src/main.o,,$(wildcard $(ATHENA_DIR)/src/*.o $(ATHENA_DIR)/src/*/*.o))

all: worker_code 

$(ATHENA_DIR):
ifdef DO_INSTALL_ATHENA
		$(DOWNLOAD_FROM_WEB) $(VERSION)
		cd $(ATHENA_DIR)/src;patch Makefile.in < ../../../makefile.patch
else
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@echo "DO_INSTALL_ATHENA is not set. Athena will not be downloaded and build."
		@echo "If you do want Athena, set DO_INSTALL_ATHENA  to 1."
		@echo "bash> export DO_INSTALL_ATHENA=1"
		@echo "csh> setenv DO_INSTALL_ATHENA  1"
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@make -s --no-print-directory -C . raise_error
endif


clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code
	$(RM) -f $(CODELIB)

codeclean: clean
	make -C $(ATHENA_DIR) clean

distclean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code
	$(RM) -Rf src

$(CODELIB): $(ATHENA_DIR) 
	cd $(ATHENA_DIR); ./configure --with-gas=hydro --enable-mpi --enable-smr
	cd $(ATHENA_DIR); ln -fs ../../../amuse_problem_$(VERSION).c src/problem.c
	cd $(ATHENA_DIR); make all CFLAGS=-g CC=mpicc LDR=mpicc


worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py AthenaInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py -i amuse.support.legacy.stopping_conditions.StoppingConditionInterface AthenaInterface -o $@

worker_code: worker_code.cc worker_code.h $(CODELIB) $(AMUSE_OBJS)
	mpicxx $(CXXFLAGS) $(SCINC) $@.cc $(AMUSE_OBJS) $(CODELIB) $(LDFLAGS) -o $@ $(SCLIBS)

.c.o: $<
	$(CC) $(CXXFLAGS) $(SCINC) -c -o $@ $< 


