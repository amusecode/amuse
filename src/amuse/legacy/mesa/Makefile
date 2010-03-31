FORTRAN ?= gfortran
export FC = $(FORTRAN)

AMUSE_DIR ?=../../../..
MESA_MAKE_DIR=./mesa_reqs
INCLUDE_DIR=./src/star/make
MESA_DIR = ./src

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
DOWNLOAD_FROM_SVN = ./build_mesa.py build

all:	worker_code

src:
ifdef DO_INSTALL_MESA
		$(DOWNLOAD_FROM_SVN)
else
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@echo "DO_INSTALL_MESA is not set. MESA will not be downloaded and build."
		@echo "If you do want MESA, set DO_INSTALL_MESA to 1."
		@echo "bash> export DO_INSTALL_MESA=1"
		@echo "csh> setenv DO_INSTALL_MESA 1"
		@echo "Note: MESA is quite large (800 MB downloaded, 2.3 GB built) and"
		@echo "requires ifort or gfortran (version > 4.3.0)."
		@echo ""
		@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@echo ""
		@make -s --no-print-directory -C . raise_error
endif

clean:
	$(RM) -f *~ *.pyc *.mod *.o worker_code worker_code.f90
	make -C $(MESA_MAKE_DIR) clean

veryclean:
	$(RM) -f *~ *.pyc *.mod *.o worker_code worker_code.f90
	make -C $(MESA_MAKE_DIR) veryclean

distclean: clean
	$(RM) -Rf ./src/

worker_code:	worker_code.f90 MESA_code.o $(STAR_WORK_OBJS)
	mpif90 $(FFLAGS) $(FCopenmp) $^ -o $@ $(LOAD_MESA_STAR)

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MESAInterface -o $@

MESA_code.o: MESA_code.f src
	make -s --no-print-directory -C $(MESA_MAKE_DIR) star FFLAGS="$(FFLAGS)"
	$(WORK_COMPILE) -I$(INCLUDE_DIR) -I$(MESA_MAKE_DIR) -o $@ $<
