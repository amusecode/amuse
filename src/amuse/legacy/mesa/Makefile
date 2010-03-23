ifdef DO_INSTALL_MESA
    TARGS = .compile
else
    TARGS = .warning
endif

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

all:	$(TARGS)

.compile:   MESA_source_files worker_code

.warning:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DO_INSTALL_MESA is not set. MESA will not be downloaded and build."
	@echo "If you do want MESA, set DO_INSTALL_MESA to any (not empty) value."
	@echo "Note: MESA is quite large (800 MB downloaded, 2.3 GB built) and"
	@echo "requires ifort or gfortran (version > 4.3.0)."
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C $(MESA_MAKE_DIR) raise_error

clean:
	/bin/rm -f *~ *.pyc *.mod *.o worker_code worker_code.f90
	make -C $(MESA_MAKE_DIR) clean

veryclean:
	/bin/rm -f *~ *.pyc *.mod *.o worker_code worker_code.f90
	make -C $(MESA_MAKE_DIR) veryclean

clean_also_mesa_checkout: clean
	rm -fr ./src/

MESA_source_files:
  ifeq ($(wildcard ./src),)
	@echo "MESA source directory not found. Will download and build MESA now..."
	./build_mesa.py build
  else
	@echo "MESA source directory was found. Continuing..."
  endif

MESA_rebuild:
	@echo "Forced rebuild of MESA. Will clean, svn-update, and re-install MESA now..."
	./build_mesa.py rebuild

worker_code:	worker_code.f90 MESA_code.o $(STAR_WORK_OBJS)
	mpif90 $(FFLAGS) $(FCopenmp) $^ -o $@ $(LOAD_MESA_STAR)

worker_code.f90: interface.py
#	@$(CODE_GENERATOR) --type=f90 interface.py MESA -o $@
	$(CODE_GENERATOR) --type=f90 interface.py MESAInterface -o $@

MESA_code.o: MESA_code.f
	make -s --no-print-directory -C $(MESA_MAKE_DIR) star FFLAGS="$(FFLAGS)"
	$(WORK_COMPILE) -I$(INCLUDE_DIR) -I$(MESA_MAKE_DIR) -o $@ $<
