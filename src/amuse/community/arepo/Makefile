# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

## examples/Noh_3d/Config.sh
## config file for 3d Noh probelm

#--------------------------------------- Basic operation mode of code
AREPOFLAGS += -DREFLECTIVE_X=2                           # in-/outflow boundary conditions in x direction
AREPOFLAGS += -DREFLECTIVE_Y=2                           # in-/outflow boundary conditions in y direction
AREPOFLAGS += -DREFLECTIVE_Z=2                           # in-/outflow boundary conditions in z direction

#--------------------------------------- Mesh motion and regularization
AREPOFLAGS += -DREGULARIZE_MESH_CM_DRIFT                 # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
AREPOFLAGS += -DREGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
AREPOFLAGS += -DREGULARIZE_MESH_FACE_ANGLE               # Use maximum face angle as roundness criterion in mesh regularization

#--------------------------------------- Time integration options
AREPOFLAGS += -DTREE_BASED_TIMESTEPS                     # non-local timestep criterion (take 'signal speed' into account)

#---------------------------------------- Single/Double Precision
AREPOFLAGS += -DDOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
AREPOFLAGS += -DINPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
AREPOFLAGS += -DOUTPUT_CENTER_OF_MASS                    # output centers of cells

#--------------------------------------- Output/Input options
AREPOFLAGS += -DHAVE_HDF5                                # needed when HDF5 I/O support is desired; should this be standard?

#--------------------------------------- Testing and Debugging options
#AREPOFLAGS += -DDEBUG                                    # enables core-dumps, should this be standard?

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
