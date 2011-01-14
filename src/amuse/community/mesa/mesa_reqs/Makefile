FORTRAN ?= gfortran
FC = $(FORTRAN)

MESA_DIR = ../src
WORK_DIR = .
STAR_TEST_DIR = $(MESA_DIR)/star/test
STAR_LIB_DIR = $(MESA_DIR)/star/make
STAR_RUN_LIB_DIR = $(MESA_DIR)/star/test/make
STAR_SRC_DIR = $(STAR_TEST_DIR)/src

include $(WORK_DIR)/makefile_header
include $(WORK_DIR)/makefile_include

LOAD_MESA = -L$(MESA_LIB_DIR) $(LOAD_MESA_STAR)
STAR_WORK_OBJS = \
   run_star_support.o run_star_extras.o \
   calibrate.o isochrone.o $(PGSTAR_OBJS) \
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

%.o: $(STAR_SRC_DIR)/%.f
	$(WORK_COMPILE) $<

clean:
	-@rm -f *.o *.mod $(WORK_DIR)/$(STAR)

remk:
	-@rm -f run.o $(WORK_DIR)/$(STAR)
