
FORTRAN ?= gfortran
AMUSE_DIR ?=../../../..
MESA_MAKE_DIR=./src/amuse
INCLUDE_DIR=./src/star/make
MESA_DIR = ./src

# FC=$(FORTRAN) # Not necessary; will be set in:
ifeq ($(wildcard ./src),)
include ./amuse/makefile_header_v1943
include ./amuse/makefile_include
else
include $(MESA_MAKE_DIR)/makefile_header_v1943
include $(MESA_MAKE_DIR)/makefile_include
endif

WORK_COMPILE = $(FC) $(FCbasic) $(FCopenmp) $(FCchecks) $(FCdebug) $(FFLAGS) \
   -I$(MESA_INCLUDE_DIR) -c $(FCfree)
STAR_WORK_OBJS = \
   $(MESA_MAKE_DIR)/run_star_support.o $(MESA_MAKE_DIR)/run_star_extras.o \
   $(MESA_MAKE_DIR)/calibrate.o $(MESA_MAKE_DIR)/isochrone.o $(MESA_MAKE_DIR)/$(NETCDF_OBJ) \
   $(MESA_MAKE_DIR)/create_zams.o $(MESA_MAKE_DIR)/sample_zams.o \
   $(MESA_MAKE_DIR)/run_star.o 

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all:   MESA_source_files worker_code

clean:
	/bin/rm -f *~ *.pyc *.mod *.o $(MESA_MAKE_DIR)/*.o $(MESA_MAKE_DIR)/star $(MESA_MAKE_DIR)/*.mod worker_code worker_code.f90 inlist

clean_also_mesa_checkout: clean
	rm -fr ./src/

MESA_source_files:
  ifeq ($(wildcard ./src),)
	echo "MESA source directory not found. Will download and build MESA now..."
	./build_mesa.py build
	cp -r ./amuse/ ./src/
  else
	echo "MESA source directory was found. Continuing..."
  endif

worker_code:	worker_code.f90 MESA_code.o $(STAR_WORK_OBJS)
	mpif90 $(FFLAGS) $(FCopenmp) $^ -o $@ $(LOAD_MESA_STAR)

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py MESA -o $@

MESA_code.o: MESA_code.f
	make -C $(MESA_MAKE_DIR) star FFLAGS="$(FFLAGS)"
	$(WORK_COMPILE) -I$(INCLUDE_DIR) -I$(MESA_MAKE_DIR) -o $@ $<
