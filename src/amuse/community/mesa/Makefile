# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

GFORTRAN_VERSION?=$(shell $(FC) -v 2>&1 | grep gcc\ version | cut -d\  -f3)
GFORTRAN_MAJOR_MINOR=$(shell echo $(GFORTRAN_VERSION) | cut -s -d'.' -f1,2)
#IFORT_VERSION?=$(shell $(FC) -v 2>&1 | cut -d\  -f2)
#IFORT_MAJOR_MINOR=$(shell echo $(IFORT_VERSION) | cut -s -d'.' -f1,2)
MPIF90 ?= mpif90
FORTRAN ?= gfortran
PYTHON ?= python
FC      = $(FORTRAN)
MPIFC+= $(FCFLAGS)
export FC MPIFC

MESA_MAKE_DIR=./mesa_reqs
INCLUDE_DIR=./src/mesa/star/make
MESA_DIR = ./src/mesa
MESA_LIB=$(MESA_DIR)/lib/libstar.a

export OPENMP_FCFLAGS=
include $(MESA_MAKE_DIR)/makefile_header
include $(MESA_MAKE_DIR)/makefile_include

WORK_COMPILE = $(FC) $(FCbasic) $(FCopenmp) $(FCchecks) $(FCdebug) $(FCFLAGS)  $(FFLAGS)\
   -I$(MESA_INCLUDE_DIR) -c $(FCfree)
STAR_WORK_OBJS = \
   $(MESA_MAKE_DIR)/run_star_support.o $(MESA_MAKE_DIR)/run_star_extras.o \
   $(MESA_MAKE_DIR)/calibrate.o $(MESA_MAKE_DIR)/isochrone.o $(MESA_MAKE_DIR)/$(PGSTAR_OBJS) \
   $(MESA_MAKE_DIR)/create_zams.o $(MESA_MAKE_DIR)/sample_zams.o \
   $(MESA_MAKE_DIR)/run_star.o 

ZAMS_DIR = $(AMUSE_DIR)/data/mesa/output/star_data/starting_models
ZAMS_DATA = $(ZAMS_DIR)/zams_z20m3.data

CODE_GENERATOR = $(AMUSE_DIR)/build.py
DOWNLOAD_FROM_WEB = $(PYTHON) ./download_svn.py
PATCH_FILES = $(PYTHON) ./patch_files.py

all:	mesa_worker $(ZAMS_DATA)

ifdef DOWNLOAD_CODES
$(MESA_DIR)/install:
	make -C . download
else
$(MESA_DIR)/install:
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@echo "DOWNLOAD_CODES is not set. MESA will not be downloaded and build."
		@echo "If you do want MESA, set DOWNLOAD_CODES to 1."
		@echo "bash> export DOWNLOAD_CODES=1"
		@echo "csh> setenv DOWNLOAD_CODES 1"
		@echo "Note: MESA is quite large (800 MB downloaded, 2.3 GB built) and"
		@echo "requires ifort (version >= 11.1) or gfortran (version >= 4.3.0)."
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@make -s --no-print-directory -C . raise_error
endif

download:
	$(RM) -Rf src
	-mkdir src
	$(DOWNLOAD_FROM_WEB)
	$(PATCH_FILES)

clean:
	$(RM) -f *~ *.pyc *.mod *.o worker_code worker_code.f90 $(MESA_LIB)
	$(RM) -f mesa_worker mesa_worker_sockets worker_code-sockets.f90
	make -C $(MESA_MAKE_DIR) clean

veryclean: clean
	cd $(MESA_DIR); ./clean

distclean: clean
	$(RM) -Rf src
	$(RM) -Rf src.*
	$(RM) -Rf .pc

$(ZAMS_DATA):
	gunzip -c $(ZAMS_DATA).gz > $(ZAMS_DATA)

mesa_worker:	worker_code.f90 interface.o $(STAR_WORK_OBJS)
	$(MPIF90) $(FCFLAGS) $(FFLAGS) $(FCopenmp) $(FS_FLAGS) $^ -o $@ $(LOAD_MESA_STAR) $(FS_LIBS)

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MESAInterface -o $@

$(MESA_LIB):  $(MESA_DIR)/install
	cp -f ./install $(MESA_DIR)
	cp -f ./test_jina.f $(MESA_DIR)/jina/test/src/test_jina.f
	cd $(MESA_DIR); ./install MPIFC="$(MPIFC) $(FCFLAGS)" OPENMP_FCFLAGS=

interface.o: interface.f $(MESA_LIB)
	make -s --no-print-directory -C $(MESA_MAKE_DIR) star FFLAGS="$(FCFLAGS)"
	$(WORK_COMPILE) -I$(INCLUDE_DIR) -I$(MESA_MAKE_DIR) -o $@ $<
