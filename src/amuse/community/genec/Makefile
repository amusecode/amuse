# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
FC      = $(MPIFC)

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

# FCFLAGS += -std=f95
# FCFLAGS += -O0

OBJS = interface.o

CODEDIR = src/GENEC/code
INI_CODEDIR = src/GENEC/tools/makeini

CODELIB = src/libgenec.a

MODULES_F90_BASE = $(CODEDIR)/evol.f90  \
                   $(CODEDIR)/const.f90  \
                   $(CODEDIR)/caramodele.f90  \
                   $(CODEDIR)/inputparam.f90  \
		   $(CODEDIR)/state.f90 \
                   $(CODEDIR)/interpolation.f90  \
                   $(CODEDIR)/nagmod.f90  \
                   $(CODEDIR)/SmallFunc.f90  \
                   $(CODEDIR)/abundmod.f90  \
                   $(CODEDIR)/diffadvmod.f90  \
                   $(CODEDIR)/equadiffmod.f90  \
                   $(CODEDIR)/geomod.f90  \
                   $(CODEDIR)/magmod.f90  \
                   $(CODEDIR)/rotmod.f90  \
                   $(CODEDIR)/strucmod.f90  \
                   $(CODEDIR)/ionisation.f90  \
                   $(CODEDIR)/convection.f90  \
                   $(CODEDIR)/timestep.f90  \
                   $(CODEDIR)/PrintAll.f90  \
                   $(CODEDIR)/omegamod.f90  \
                   $(CODEDIR)/EOS.f90  \
                   $(CODEDIR)/nablas.f90  \
                   $(CODEDIR)/energy.f90  \
                   $(CODEDIR)/opacity.f90  \
                   $(CODEDIR)/envelope.f90  \
                   $(CODEDIR)/LayersShift.f90  \
                   $(CODEDIR)/advection.f90  \
                   $(CODEDIR)/diffusion.f90  \
                   $(CODEDIR)/chemicals.f90  \
                   $(CODEDIR)/bintidemod.f90 \
                   $(CODEDIR)/winds.f90 \
		   $(INI_CODEDIR)/inichemmod.f90 \
                   $(INI_CODEDIR)/interpolmod.f90 \
                   $(INI_CODEDIR)/modinimod.f90

MODULES_F90 = $(MODULES_F90_BASE) \
              $(CODEDIR)/PGPlotModule_Null.f90  \
              $(CODEDIR)/henyey_solver.f90  \
              $(CODEDIR)/WriteSaveClose.f90 \
	      $(CODEDIR)/genec.f90
# SOURCES_F90 = $(CODEDIR)/main.f90


MODULES_F90_ALL = $(MODULES_F90) \
		  $(MODULES_F90_MAKEINI) \
		  $(CODEDIR)/amuse.f90

MODULES=$(MODULES_F90_ALL:%.f90=%.o)
# SOURCES=$(SOURCES_F90:%.f90=%.o)

OBJECTS=$(MODULES) # $(SOURCES)

all: genec_worker 

clean:
	$(RM) -rf __pycache__
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ genec_worker worker_code.f90
	$(RM) $(CODELIB)
	$(RM) $(CODEDIR)/*.o
	$(RM) $(INI_CODEDIR)/*.o
	$(RM) *.mod
	# make -C $(CODEDIR) -f make_amuse clean

distclean: clean
	make -C src distclean

$(CODELIB): $(OBJECTS)
	@echo "==========LINKING========="
	@echo "Fortran options: $(FFLAGS) $(LDFLAGS) --> $(CODELIB)"
	ar rc $(CODELIB) $(OBJECTS)
	# @$(COMPILER) -shared $(FFLAGS) $(OBJECTS) $(LDFLAGS) -o $(CODELIB)


worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py GenecInterface -o $@

genec_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $< $(OBJS) $(CODELIB) $(FS_LIBS) -o $@

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<
