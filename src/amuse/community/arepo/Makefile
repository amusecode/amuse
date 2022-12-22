# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

AREPOFLAGS += -DSELFGRAVITY                              # gravitational intraction between simulation particles/cells
AREPOFLAGS += -DHIERARCHICAL_GRAVITY                     # use hierarchical splitting of the time integration of the gravity
AREPOFLAGS += -DCELL_CENTER_GRAVITY                      # uses geometric centers to calculate gravity of cells, only possible with HIERARCHICAL_GRAVITY
AREPOFLAGS += -DALLOW_DIRECT_SUMMATION                   # Performed direct summation instead of tree-based gravity if number of active particles < DIRECT_SUMMATION_THRESHOLD (= 3000 unless specified differently here)
AREPOFLAGS += -DDIRECT_SUMMATION_THRESHOLD=500           # Overrides maximum number of active particles for which direct summation is performed instead of tree based calculation
AREPOFLAGS += -DGRAVITY_NOT_PERIODIC                     # gravity is not treated periodically
AREPOFLAGS += -DNSOFTTYPES=2                             # Number of different softening values to which particle types can be mapped.
AREPOFLAGS += -DMULTIPLE_NODE_SOFTENING                  # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
AREPOFLAGS += -DINDIVIDUAL_GRAVITY_SOFTENING=32          # bitmask with particle types where the softenig type should be chosen with that of parttype 1 as a reference type
AREPOFLAGS += -DADAPTIVE_HYDRO_SOFTENING                 # Adaptive softening of gas cells depending on their size
AREPOFLAGS += -DTREE_BASED_TIMESTEPS                     # non-local timestep criterion (take 'signal speed' into account)
AREPOFLAGS += -DDOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
AREPOFLAGS += -DNGB_TREE_DOUBLEPRECISION                 # if this is enabled, double precision is used for the neighbor node extension
AREPOFLAGS += -DPROCESS_TIMES_OF_OUTPUTLIST              # goes through times of output list prior to starting the simulaiton to ensure that outputs are written as close to the desired time as possible (as opposed to at next possible time if this flag is not active)
AREPOFLAGS += -DHAVE_HDF5                                # needed when HDF5 I/O support is desired (recommended)
# AREPOFLAGS += -DDEBUG                                    # enables core-dumps

HDF5_FLAGS += -DH5_USE_16_API

MPICXX   ?= mpicxx
CXX = $(MPICXX)
CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) $(GSL_FLAGS) $(AREPOFLAGS) $(HDF5_FLAGS)
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libarepo.a

all: arepo_worker 

clean:
	$(RM) -rf __pycache__
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ arepo_worker worker_code.cc
	make -C src clean

distclean: clean
	make -C src distclean

$(CODELIB): .FORCE
	export AREPOFLAGS
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py ArepoInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py ArepoInterface -o $@

arepo_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(GSL_FLAGS) $(HDF5_FLAGS) $< $(OBJS) $(CODELIB) -o $@ $(GMP_LIBS) $(GSL_LIBS) $(HDF5_LIBS)

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) $(CODELIB) -c -o $@ $< 

.FORCE:
