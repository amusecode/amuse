# standard amuse configuration include
# config.mk will be made after ./configure has run
#ifeq ($(origin AMUSE_DIR), undefined)
AMUSE_DIR := $(shell amusifier --get-amuse-dir)
#endif
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

use_phantom_grape_x86 = no  # speedup from ~1M particles or so?

FDPS_PATH = src/FDPS/src
SUI_PATH = src/sui

INCLUDE += -I $(FDPS_PATH) -I $(SUI_PATH)

MPIFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

CXXFLAGS = $(INCLUDE) -fPIC -O3 -Wall -std=c++20 -fopenmp


CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
CXXFLAGS += -DENABLE_HYDRO_INTERACT
# CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DENABLE_VARIABLE_SMOOTHING_LENGTH
# CXXFLAGS += -DUSE_ENTROPY
# CXXFLAGS += -DUSE_BALSARA_SWITCH
# CXXFLAGS += -DUSE_PRESCR_OF_THOMAS_COUCHMAN_1992
CXXFLAGS += -DISOTHERMAL_EOS

ifeq ($(use_phantom_grape_x86),yes)
  PG_ROOT = $(FDPS_PATH)/phantom_grape_x86/G5/newton/libpg5
  CXXFLAGS += -I $(PG_ROOT)
  CXXFLAGS += -DENABLE_PHANTOM_GRAPE_X86
  CXXFLAGS += -L$(PG_ROOT) -lpg5
  PG_BUILD = make -C $(PG_ROOT) distclean libpg5.a
  PG_CLEAN = make -C $(PG_ROOT) distclean
else
  PG_BUILD =
  PG_CLEAN = 
endif

LDFLAGS += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o


all: sui_worker $(CODELIB)

$(CODELIB):
	$(PG_BUILD)
	#amusi	$(MPICXX) $(CXXFLAGS)
	#$(PG_BUILD)
	#$(AR) $@ 

clean:
	$(RM) -rf __pycache__
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h
	$(RM) *~ sui_worker worker_code.cc
	make -C src clean

distclean: clean
	make -C src distclean

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py SuiInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H  -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py SuiInterface -o $@

sui_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS) interface.h
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(LDFLAGS) $< $(OBJS) -o $@

interface.o: interface.cc
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) -c -o $@ $<
