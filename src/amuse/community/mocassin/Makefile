# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPIF90 ?= mpif90
PYTHON ?= python
FC      = $(MPIF90)


OBJS = interface.o

CODE_GENERATOR = $(AMUSE_DIR)/build.py

VERSION = 2.02.69
CODEDIR = src/mocassin.$(VERSION)
CODELIB = $(CODEDIR)/libmocassin.a

SOURCEDIR = $(CODEDIR)/source
SOURCES  = \
	$(SOURCEDIR)/infnan.f90 $(SOURCEDIR)/constants_mod.f90 $(SOURCEDIR)/vector_mod.f90 \
	$(SOURCEDIR)/common_mod.f90 $(SOURCEDIR)/interpolation_mod.f90 $(SOURCEDIR)/set_input_mod.f90 \
	$(SOURCEDIR)/hydro_mod.f90 $(SOURCEDIR)/ph_mod.f90 $(SOURCEDIR)/composition_mod.f90 \
	$(SOURCEDIR)/continuum_mod.f90 $(SOURCEDIR)/ionization_mod.f90 $(SOURCEDIR)/pathIntegration_mod.f90 \
	$(SOURCEDIR)/grid_mod.f90 $(SOURCEDIR)/dust_mod.f90 $(SOURCEDIR)/emission_mod.f90 $(SOURCEDIR)/photon_mod.f90  \
	$(SOURCEDIR)/update_mod.f90 $(SOURCEDIR)/output_mod.f90 $(SOURCEDIR)/iteration_mod.f90
	
OBJECTS = $(SOURCES:.f90=.o)


DOWNLOAD_FROM_WEB = $(PYTHON) ./downloadhttp.py
PATCH_FILES = $(PYTHON) ./patch_files.py


FFLAGS   += -Wall -g -I$(CODEDIR)/source
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

FSOPTS ?= -L$(AMUSE_DIR)/lib/forsockets -lforsockets -I$(AMUSE_DIR)/lib/forsockets

all: mocassin_worker 

ifdef DOWNLOAD_CODES
$(CODEDIR)/makefile:
	make -C . download
else
$(CODEDIR)/makefile:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. Mocassin will not be downloaded and build."
	@echo "If you do want Mocassin, set DOWNLOAD_CODES to 1."
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
	$(RM) -f *.mod *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code worker_code.f90
	$(RM) mocassin_worker
	-make -C $(CODEDIR) clean
	
distclean: clean
	$(RM) -Rf src
	$(RM) -Rf .pc
	$(RM) -Rf src.*

$(CODELIB): $(CODEDIR)/makefile
	make -C $(CODEDIR) clean
	make -C $(CODEDIR) all

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MocassinInterface -o $@

mocassin_worker:  worker_code.f90 $(CODELIB) $(OBJS) 
	$(MPIF90) $(FFLAGS) worker_code.f90 $(OBJS) $(OBJECTS) -o $@
	
mocassin_worker_sockets:  worker_code.f90 $(CODELIB) $(OBJS)  
	$(MPIF90) $(FFLAGS) worker_code.f90 $(OBJS) $(OBJECTS) $(FSOPTS) -o $@

%.mod %.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<
