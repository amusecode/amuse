FORTRAN ?= gfortran
FC = $(FORTRAN)

MESA_DIR = ../src
WORK_DIR = .
WORK_SRC_DIR = $(WORK_DIR)/src
STAR_TEST_DIR = $(MESA_DIR)/star/test
STAR_LIB_DIR = $(MESA_DIR)/star/make
STAR_RUN_LIB_DIR = $(MESA_DIR)/star/test/make
STAR_SRC_DIR = $(STAR_TEST_DIR)/src

include $(WORK_DIR)/makefile_header_v1943
include $(WORK_DIR)/makefile_include

LOAD_MESA = -L$(MESA_LIB_DIR) $(LOAD_MESA_STAR)
STAR_WORK_OBJS = \
   run_star_support.o run_star_extras.o \
   calibrate.o isochrone.o $(NETCDF_OBJ) \
   create_zams.o sample_zams.o \
   run_star.o 

#################################################################

ifndef STAR
STAR = star
endif

OBJS = $(STAR_WORK_OBJS) run.o

$(STAR) : lib_mesa_star $(OBJS)
	$(FC) $(FCopenmp) $(FFLAGS) -o $(WORK_DIR)/$(STAR) $(OBJS) $(LOAD_MESA_STAR)
	
#################################################################

# change this as necessary.  see utils/makefile_header for definitions.
WORK_COMPILE = $(FC) $(FCbasic) $(FCopenmp) $(FCchecks) $(FCdebug) $(FFLAGS) \
   -I$(MESA_INCLUDE_DIR) -c $(FCfree)

lib_mesa_star:
	make -C $(STAR_LIB_DIR) libstar.a FFLAGS="$(FFLAGS)"
	cp $(STAR_LIB_DIR)/star_lib.mod $(MESA_DIR)/include
	cp $(STAR_LIB_DIR)/star_def.mod $(MESA_DIR)/include
	cp $(STAR_LIB_DIR)/libstar.a $(MESA_DIR)/lib
	ranlib $(MESA_DIR)/lib/libstar.a

run.o: $(WORK_SRC_DIR)/run.f
	$(WORK_COMPILE) $<

run_star_extras.o: $(WORK_SRC_DIR)/run_star_extras.f
	$(WORK_COMPILE) $<

create_zams.o: $(WORK_SRC_DIR)/create_zams.f
	$(WORK_COMPILE) $<

%.o: $(STAR_SRC_DIR)/%.f
	$(WORK_COMPILE) $<

clean:
	-@rm -f *.o *.mod $(WORK_DIR)/$(STAR)

veryclean: clean
	make -C $(STAR_LIB_DIR) clean
	make -C $(STAR_RUN_LIB_DIR) clean

remk:
	-@rm -f run.o $(WORK_DIR)/$(STAR)
