MPICXX ?= mpicxx
MPICC ?= mpicc
# TODO: Determine whether this is needed as included in arepo/Makefile.
ifeq ($(origin AMUSE_DIR), undefined)
	AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

CC       =  $(MPICC)  # sets the C-compiler

# GSL_FLAGS  = -I/opt/Homebrew/include  # TODO: Need to make GSL_FLAGS generalisable.
CFLAGS   += -Wall -g $(GSL_FLAGS)
CXXFLAGS += $(CFLAGS) 

HDF5_FLAGS += -DH5_USE_16_API

LDFLAGS  += -lm $(MUSE_LD_FLAGS) $(HDF5_FLAGS)
LIBS += $(HDF5_LIBS)

CODELIB = libarepo.a

SUBDIR = add_backgroundgrid cooling debug_md5 domain fof gitversion gravity \
	hydro init io mesh mesh/voronoi mpi_utils ngbtree star_formation subfind \
	time_integration utils
SRCS = $(foreach fd, $(SUBDIR), $(wildcard $(fd)/*.c))
SRCS += main/allvars.c

CODEOBJS = $(SRCS:c=o)
CODEOBJS += main/allvars.o main/run.o main/main.o

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


CXXFLAGS += $(AREPOFLAGS)

# Add includes
CXXFLAGS += $(HDF5_FLAGS)

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB)

clean:
	$(RM) -f *.o *.a
	$(RM) $(SRCS:c=o)

distclean: clean

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< $(LIBS) 

.c.o: $<
	$(MPICC) $(CXXFLAGS) -c -o $@ $< $(LIBS) 
