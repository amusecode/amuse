MPIF90 ?= mpif90
FORTRAN ?= gfortran
FC      = $(FORTRAN)
export FC

AMUSE_DIR ?=../../../..
MESA_MAKE_DIR=./mesa_reqs
INCLUDE_DIR=./src/star/make
MESA_DIR = ./src
MESA_LIB=$(MESA_DIR)/lib/libstar.a

FSLIBS ?= -L$(AMUSE_DIR)/lib/forsockets -lforsockets
FSINC ?= -I$(AMUSE_DIR)/lib/forsockets

include $(MESA_MAKE_DIR)/makefile_header
include $(MESA_MAKE_DIR)/makefile_include

WORK_COMPILE = $(FC) $(FCbasic) $(FCopenmp) $(FCchecks) $(FCdebug) $(FFLAGS) \
   -I$(MESA_INCLUDE_DIR) -c $(FCfree)
STAR_WORK_OBJS = \
   $(MESA_MAKE_DIR)/run_star_support.o $(MESA_MAKE_DIR)/run_star_extras.o \
   $(MESA_MAKE_DIR)/calibrate.o $(MESA_MAKE_DIR)/isochrone.o $(MESA_MAKE_DIR)/$(PGSTAR_OBJS) \
   $(MESA_MAKE_DIR)/create_zams.o $(MESA_MAKE_DIR)/sample_zams.o \
   $(MESA_MAKE_DIR)/run_star.o 

CODE_GENERATOR = $(AMUSE_DIR)/build.py
DOWNLOAD_FROM_SVN = ./build_mesa.py download

all:	mesa_worker


ifdef DOWNLOAD_CODES
src: download
else
src:
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
	$(DOWNLOAD_FROM_SVN)

clean:
	$(RM) -f *~ *.pyc *.mod *.o worker_code worker_code.f90 $(MESA_LIB)
	$(RM) -f mesa_worker mesa_worker_sockets worker_code-sockets.f90
	make -C $(MESA_MAKE_DIR) clean

veryclean: clean
	cd $(MESA_DIR); ./clean

distclean: clean
	$(RM) -Rf ./src/

mesa_worker:	worker_code.f90 MESA_code.o $(STAR_WORK_OBJS)
	$(MPIF90) $(FFLAGS) $(FCopenmp) $^ -o $@ $(LOAD_MESA_STAR)

mesa_worker_sockets: worker_code-sockets.f90 MESA_code.o $(STAR_WORK_OBJS)
	$(MPIF90) $(FFLAGS) $(FCopenmp) $(FSINC) $^ -o $@ $(LOAD_MESA_STAR) $(FSLIBS)

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MESAInterface -o $@

worker_code-sockets.f90: interface.py
	$(CODE_GENERATOR) --type=f90 --mode=sockets interface.py MESAInterface -o $@

$(MESA_LIB): src
	cd $(MESA_DIR); ./install

MESA_code.o: MESA_code.f $(MESA_LIB)
	make -s --no-print-directory -C $(MESA_MAKE_DIR) star FFLAGS="$(FFLAGS)"
	$(WORK_COMPILE) -I$(INCLUDE_DIR) -I$(MESA_MAKE_DIR) -o $@ $<
