
ATHENA_DIR = src/athena3.1

CC = mpicc
CFLAGS   += -g -I$(ATHENA_DIR)/src 
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -ldl -lm

AMUSE_OBJS = interface.o

CODELIB = libathena3.1.a

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py
DOWNLOAD_FROM_WEB = ./downloadhttp.py

ATHENA_OBJS := $(subst $(ATHENA_DIR)/src/main.o,,$(wildcard $(ATHENA_DIR)/src/*.o))

all: worker_code 

src:
ifdef DO_INSTALL_ATHENA
		$(DOWNLOAD_FROM_WEB)
else
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@echo "DO_INSTALL_ATHENA is not set. Athena will not be downloaded and build."
		@echo "If you do want SeBa, set DO_INSTALL_ATHENA  to 1."
		@echo "bash> export DO_INSTALL_ATHENA = 1"
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
	make -C src clean
	
distclean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code
	$(RM) -Rf src

$(CODELIB): src
	@echo $(ATHENA_OBJS)
	cd $(ATHENA_DIR); ./configure --with-gas=hydro --enable-mpi
	cd $(ATHENA_DIR); ln -fs ../../../amuse_problem.c src/problem.c
	cd $(ATHENA_DIR); make all CFLAGS=-g CC=mpicc LDR=mpicc
	ar rv $(CODELIB) $(ATHENA_OBJS)
	ranlib $(CODELIB)
     
worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py AthenaInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py AthenaInterface -o $@

worker_code: worker_code.cc worker_code.h $(CODELIB) $(AMUSE_OBJS)
	mpicxx $(CXXFLAGS) $@.cc $(AMUSE_OBJS) $(CODELIB) $(LDFLAGS) -o $@

.c.o: $<
	$(CC) $(CXXFLAGS) -c -o $@ $< 
	
	
