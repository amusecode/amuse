MESA_DIR = ../src
WORK_DIR = .
WORK_SRC_DIR = $(WORK_DIR)/src
STAR_TEST_DIR = $(MESA_DIR)/star/test
STAR_SRC_DIR = $(STAR_TEST_DIR)/src

#include $(MESA_DIR)/utils/makefile_header
#include $(MESA_DIR)/star/test/make/makefile_include
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

$(STAR) : $(OBJS)
	$(FC) $(FCopenmp) $(FFLAGS) -o $(WORK_DIR)/$(STAR) $(OBJS) $(LOAD_MESA_STAR)
	
#################################################################

# change this as necessary.  see utils/makefile_header for definitions.
WORK_COMPILE = $(FC) $(FCbasic) $(FCopenmp) $(FCchecks) $(FCdebug) $(FFLAGS) \
   -I$(MESA_INCLUDE_DIR) -c $(FCfree)

run.o: $(WORK_SRC_DIR)/run.f
	$(WORK_COMPILE) $<

run_star_extras.o: $(WORK_SRC_DIR)/run_star_extras.f
	$(WORK_COMPILE) $<

%.o: $(STAR_SRC_DIR)/%.f
	$(WORK_COMPILE) $<

clean:
	-@rm -f *.o *.mod $(WORK_DIR)/$(STAR)

remk:
	-@rm -f run.o $(WORK_DIR)/$(STAR)
