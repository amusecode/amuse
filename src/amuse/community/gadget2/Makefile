MPICXX ?= mpicxx

SRCDIR		= src
include $(SRCDIR)/makefile_options
PERIODIC_OPT := "$(OPT) -DPERIODIC -DNOGRAVITY"

CXXFLAGS += -Wall -g -DTOOLBOX $(OPT)

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

SCLIBS ?= -L../../../../lib/stopcond -lstopcond
SCINC ?= -I../../../../lib/stopcond

CODELIB = $(SRCDIR)/libgadget.a

ALLVARS = $(SRCDIR)/allvars.o

LIBFILES = $(SRCDIR)/run.c  $(SRCDIR)/predict.c  $(SRCDIR)/begrun.c \
	$(SRCDIR)/endrun.c  $(SRCDIR)/global.c  $(SRCDIR)/timestep.c  \
	$(SRCDIR)/init.c  $(SRCDIR)/restart.c  $(SRCDIR)/io.c  \
	$(SRCDIR)/accel.c  $(SRCDIR)/read_ic.c  $(SRCDIR)/ngb.c  \
	$(SRCDIR)/system.c  $(SRCDIR)/allocate.c  $(SRCDIR)/density.c  \
	$(SRCDIR)/gravtree.c  $(SRCDIR)/hydra.c  $(SRCDIR)/driftfac.c  \
	$(SRCDIR)/domain.c  $(SRCDIR)/allvars.c  $(SRCDIR)/potential.c  \
	$(SRCDIR)/forcetree.c  $(SRCDIR)/peano.c  $(SRCDIR)/gravtree_forcetest.c  \
	$(SRCDIR)/pm_periodic.c  $(SRCDIR)/pm_nonperiodic.c  $(SRCDIR)/longrange.c  \
	$(SRCDIR)/density_at_pos.c

INTERFACE_OBJ = gadget_code.o
INTERFACE_SRC  = $(INTERFACE_OBJ:.o=.cc)

all:	gadget2_worker

$(CODELIB): $(LIBFILES) $(SRCDIR)/Makefile $(SRCDIR)/makefile_options
	make -C $(SRCDIR) all

clean:
	rm -f *.o *.pyc worker_code.cc worker_code.h worker_code worker_code_periodic
	rm -f gadget2_worker gadget2_worker_periodic
	make -C $(SRCDIR) clean

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c $^ Gadget2Interface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H $^ -i amuse.support.codes.stopping_conditions.StoppingConditionInterface Gadget2Interface -o $@

gadget2_worker: worker_code.cc worker_code.h $(CODELIB) $(INTERFACE_OBJ) $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) $(SCINC) $< $(INTERFACE_OBJ) -o $@  $(CODELIB) $(SCLIBS) -lgsl -lgslcblas $(ALLVARS)
	
gadget2_worker_periodic: worker_code.cc worker_code.h $(INTERFACE_SRC) $(SRCDIR)/makefile_options
	rm -f $(INTERFACE_OBJ)
	make $(INTERFACE_OBJ) OPT=$(PERIODIC_OPT)
	make -C $(SRCDIR) clean 
	make -C $(SRCDIR) all OPT=$(PERIODIC_OPT)
	$(MPICXX) $(CXXFLAGS) -DPERIODIC $(SCINC) $< $(INTERFACE_OBJ) -o $@  $(CODELIB) $(SCLIBS) -lgsl -lgslcblas $(ALLVARS)
	rm -f $(INTERFACE_OBJ)
	make -C $(SRCDIR) clean 

gadget_code.o: gadget_code.cc $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) $(SCINC) -c -o $@ gadget_code.cc $(SCLIBS)

