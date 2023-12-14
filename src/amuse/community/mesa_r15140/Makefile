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
MPIFC ?= mpif90
MPIF90 ?= mpif90
MPICC ?= mpicc
FORTRAN ?= gfortran
PYTHON ?= python
FC      ?= $(FORTRAN)
export CC		?= gcc
MPIFC+= $(FCFLAGS)
MPICC+= $(CFLAGS)

LAPACK_LIBS ?= -llapack
BLAS_LIBS ?= -lblas

export FC MPIFC MPICC LAPACK_LIBS BLAS_LIBS


SRC_DIR := $(realpath ./src)

MESA_MAKE_DIR=./mesa_reqs
INCLUDE_DIR=./src/mesa/star/make
export MESA_DIR = $(realpath ./src/mesa-r15140)
MESA_LIB=$(MESA_DIR)/lib/libstar.a
WORK_SRC_DIR=./

export OPENMP_FCFLAGS=$(OPENMP_FCFLAG)

DOWNLOAD_FROM_WEB = $(PYTHON) ./download.py
BUILD_DEP = $(PYTHON) ./build-dep.py
PATCHES_DIR=./patches

include $(PATCHES_DIR)/makefile_header_non_mesasdk

all:	mesa_worker 

$(MESA_DIR)/install:
	make -C . download

download:
	$(RM) -Rf src
	-mkdir src
	$(DOWNLOAD_FROM_WEB)
	make -C . patch
	make -C . skip_tests
	make -C . build_hdf5
	make -C . build_ndiff
	make -C . build_makedepf90
	make -C . build_crmath
	make -C . build_crlibm
	make -C . build_fpx3_deps
	make -C . build_fpx3
	make -C . build_lapack95
	make -C . install_deps

build_ndiff:
	cd $(MESA_DIR)/utils;\
	rm -rf ndiff-2.00;\
	tar xvf ndiff-2.00.tar.gz;\
	cd ndiff-2.00;\
	./configure;\
	make;\

build_makedepf90:
	cd $(MESA_DIR)/utils;\
	rm -rf makedepf90-2.8.8;\
	tar xvf makedepf90-2.8.8.tar.gz;\
	cd 	makedepf90-2.8.8;\
	./configure;\
	make;\

build_crmath:
	cd $(SRC_DIR)/crmath-1.2;\
	make;\

build_crlibm:
	cd $(SRC_DIR)/crlibm-1.0beta4;\
	./configure --enable-static=yes;\
	make;\
	make install;\
	ar cru libcrlibm.a *.o scs_lib/*.o

build_fpx3:
	mv $(SRC_DIR)/fpx3 $(SRC_DIR)/fpx3_folder;\
	cp $(SRC_DIR)/fpx3_folder/fpx3 $(SRC_DIR)/fpx3

build_lapack95:
	cd $(SRC_DIR)/LAPACK95/SRC;\
	mkdir ../lapack95_modules;\
	make clean;\
	make single_double_complex_dcomplex

build_fpx3_deps:
	chmod u+x $(SRC_DIR)/fpx3_deps

build_hdf5:
	cd $(SRC_DIR)/hdf5-1.12.0;\
	./configure --enable-fortran --enable-fortran2003 --enable-parallel --prefix=$(MESA_DIR) --enable-shared=no;\
	make;\
	make install

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
	cd $(MESA_DIR);\
	env PATH=$(MESA_DIR)/../:$(PATH) ./install
 

patch:
	mkdir -p $(MESA_DIR)/lib $(MESA_DIR)/include
	$(PYTHON) ./patch_files.py
	cp $(PATCHES_DIR)/makefile_header_non_mesasdk $(MESA_DIR)/utils/makefile_header

install_deps:
	cp ./src/crmath-1.2/libcrmath.a $(MESA_DIR)/lib/.
	cp ./src/crmath-1.2/*.mod $(MESA_DIR)/include/.
	cp ./src/crmath-1.2/*.smod $(MESA_DIR)/include/.
	cp ./src/crlibm-1.0beta4/libcrlibm.a $(MESA_DIR)/lib/.
	cp ./src/LAPACK95/lapack95_modules/*.mod $(MESA_DIR)/include/.
	cp ./src/LAPACK95/lapack95.a $(MESA_DIR)/lib/liblapack95.a
	cp $(MESA_DIR)/utils/makedepf90-2.8.8/makedepf90 $(SRC_DIR)/.
	cp $(MESA_DIR)/utils/ndiff-2.00/ndiff $(SRC_DIR)/.

interface.o: interface.f90 $(MESA_LIB)
	$(WORK_COMPILE) $<

skip_tests:
#This one should always be skipped
	touch $(MESA_DIR)/adipls/skip_test
# We cant guarentee bit-for-bit as we dont control the compiler used
# So everyone has to skip the build tests
	touch $(MESA_DIR)/skip_test
