# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk

GFORTRAN_VERSION?=$(shell $(FC) -v 2>&1 | grep gcc\ version | cut -d\  -f3)
GFORTRAN_MAJOR_MINOR=$(shell echo $(GFORTRAN_VERSION) | cut -s -d'.' -f1,2)
#IFORT_VERSION?=$(shell $(FC) -v 2>&1 | cut -d\  -f2)
#IFORT_MAJOR_MINOR=$(shell echo $(IFORT_VERSION) | cut -s -d'.' -f1,2)
MPIF90 ?= mpif90
MPICC ?= mpicc
FORTRAN ?= gfortran
PYTHON ?= python
FC      ?= $(FORTRAN)
export CC		?= gcc
MPIFC+= $(FCFLAGS)
MPICC+= $(CFLAGS)


export FC MPIFC MPICC


MESA_MAKE_DIR=./mesa_reqs
INCLUDE_DIR=./src/mesa/star/make
export MESA_DIR = $(realpath ./src/mesa-r15140)
MESA_LIB=$(MESA_DIR)/lib/libstar.a
WORK_SRC_DIR=./

export OPENMP_FCFLAGS=-fopenmp

CODE_GENERATOR ?= $(AMUSE_DIR)/build.py
DOWNLOAD_FROM_WEB = $(PYTHON) ./download.py
BUILD_DEP = $(PYTHON) ./build-dep.py
PATCHES_DIR=./patches

include $(PATCHES_DIR)/makefile_header_non_mesasdk


all:	mesa_worker 

ifdef DOWNLOAD_CODES
$(MESA_DIR)/install:
	make -C . download
else
$(MESA_DIR)/install:
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@echo "DOWNLOAD_CODES is not set. MESA will not be downloaded and built."
		@echo "If you do want MESA, set DOWNLOAD_CODES to 1."
		@echo "bash> export DOWNLOAD_CODES=1"
		@echo "csh> setenv DOWNLOAD_CODES 1"
		@echo "Note: MESA is quite large (2.5 GB downloaded, 19 GB built) and"
		@echo "requires gfortran (version >= 9.0.0)."
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@make -s --no-print-directory -C . raise_error
endif

download:
	$(RM) -Rf src
	-mkdir src
	$(DOWNLOAD_FROM_WEB)
	make -C . build_ndiff
	make -C . build_makedepf90
	make -C . patch

build_ndiff:
	cd $(MESA_DIR)/utils;\
	rm -rf ndiff-2.00;\
	tar xvf ndiff-2.00.tar.gz;\
	cd ndiff-2.00;\
	./configure;\
	make

build_makedepf90:
	cd $(MESA_DIR)/utils;\
	rm -rf makedepf90-2.8.8;\
	tar xvf makedepf90-2.8.8.tar.gz;\
	cd 	makedepf90-2.8.8;\
	./configure;\
	make

clean:
	$(RM) -f *~ *.pyc *.mod *.smod *.o worker_code worker_code.f90 $(MESA_LIB)
	$(RM) -f mesa_worker mesa_worker_sockets worker_code-sockets.f90
	make -C $(MESA_MAKE_DIR) clean

veryclean: clean
	cd $(MESA_DIR); ./clean

distclean: clean
	$(RM) -Rf src
	$(RM) -Rf src.*
	$(RM) -Rf .pc

mesa_worker:	worker_code.f90  mesa_interface.o interface.o run_star_extras.o
	$(MPIF90) $(FCFLAGS) $(FFLAGS) $(LDFLAGS) $(FCopenmp) $(FS_FLAGS) $^ -o $@ \
		$(LOAD_EXTRAS1) $(OTHER_INCLUDES) $(LOAD_MESA_STAR) $(LOAD_EXTRAS) $(FS_LIBS)  $(LIBS)

mesa_interface.o: mesa_interface.f90 $(MESA_LIB) run_star_extras.o
	$(WORK_COMPILE) -I$(MESA_DIR)/star/make -I./ $<


run_star_extras.o: run_star_extras.f90 $(MESA_LIB)
	$(WORK_COMPILE) $<


worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MESAInterface -o $@


$(MESA_LIB):  $(MESA_DIR)/install
	mkdir -p $(MESA_DIR)/lib $(MESA_DIR)/include
	cd $(MESA_DIR); ./install
 

patch:
	mkdir -p $(MESA_DIR)/lib $(MESA_DIR)/include
	$(PYTHON) ./patch_files.py
	cp $(MESA_DIR)/utils/makefile_header_non_mesasdk $(MESA_DIR)/utils/makefile_header
	cp ./src/crmath-1.2/libcrmath.a $(MESA_DIR)/lib/.
	cp ./src/crmath-1.2/*.mod $(MESA_DIR)/include/.
	cp ./src/crmath-1.2/*.smod $(MESA_DIR)/include/.
	cp ./src/crlibm-1.0beta4/libcrlibm.a $(MESA_DIR)/lib/.


interface.o: interface.f90 $(MESA_LIB)
	$(WORK_COMPILE) $<
