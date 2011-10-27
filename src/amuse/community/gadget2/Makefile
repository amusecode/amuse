# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx

SRCDIR		= src
include $(SRCDIR)/makefile_options


BUILDDIR = build
PERIODIC_BUILDDIR = build_periodic

PERIODIC_OPT := $(OPT) -DPERIODIC -DNOGRAVITY

CXXFLAGS += -Wall -g -DTOOLBOX $(OPT)

CODE_GENERATOR = $(AMUSE_DIR)/build.py

SCLIBS ?= -L../../../../lib/stopcond -lstopcond
SCINC ?= -I../../../../lib/stopcond

GSL_FLAGS ?= $(shell gsl-config --cflags)
GSL_LIBS ?= $(shell gsl-config --libs)
GSL_INCL = $(GSL_FLAGS)


CODELIB = $(BUILDDIR)/libgadget.a

ALLVARS = $(BUILDDIR)/allvars.o

LIBFILES = $(SRCDIR)/run.c  $(SRCDIR)/predict.c  $(SRCDIR)/begrun.c \
	$(SRCDIR)/endrun.c  $(SRCDIR)/global.c  \
	$(SRCDIR)/init.c  $(SRCDIR)/restart.c  $(SRCDIR)/io.c  \
	$(SRCDIR)/accel.c  $(SRCDIR)/read_ic.c  $(SRCDIR)/ngb.c  \
	$(SRCDIR)/system.c  $(SRCDIR)/allocate.c  $(SRCDIR)/density.c  \
	$(SRCDIR)/gravtree.c  $(SRCDIR)/hydra.c  $(SRCDIR)/driftfac.c  \
	$(SRCDIR)/domain.c  $(SRCDIR)/allvars.c  $(SRCDIR)/potential.c  \
	$(SRCDIR)/forcetree.c  $(SRCDIR)/peano.c  $(SRCDIR)/gravtree_forcetest.c  \
	$(SRCDIR)/pm_periodic.c  $(SRCDIR)/pm_nonperiodic.c  $(SRCDIR)/longrange.c  \
	$(SRCDIR)/density_at_pos.c \
	$(SRCDIR)/timestep_limit.c  $(SRCDIR)/sigvel.c
#  ^^^ Durier & Dalla Vecchia timestepping fix ^^^

#	$(SRCDIR)/timestep.c
#  ^^^ original timestepping ^^^

INTERFACE_OBJ = gadget_code.o
INTERFACE_SRC  = $(INTERFACE_OBJ:.o=.cc)

all:	gadget2_worker


clean:
	rm -f *.o *.po *.pyc worker_code.cc worker_code.h worker_code worker_code_periodic
	rm -f gadget2_worker gadget2_worker_periodic gadget2_worker_sockets worker_code-sockets.cc
	rm -f gadget2_worker_periodic_sockets
	make -C $(SRCDIR) clean
	rm -Rf $(PERIODIC_BUILDDIR)
	rm -Rf $(BUILDDIR)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c $^ Gadget2Interface -o $@ 

worker_code-sockets.cc: interface.py
	$(CODE_GENERATOR) --type=c --mode=sockets $^ Gadget2Interface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H $^ -i amuse.support.codes.stopping_conditions.StoppingConditionInterface Gadget2Interface -o $@


$(BUILDDIR)/Makefile: $(SRCDIR)/Makefile 
	-mkdir $(BUILDDIR) 
	cp $(SRCDIR)/Makefile $(BUILDDIR)/Makefile
	
$(BUILDDIR)/makefile_options: $(SRCDIR)/makefile_options 
	-mkdir $(BUILDDIR) 
	cp $(SRCDIR)/makefile_options $(BUILDDIR)/makefile_options
	
	
$(CODELIB): $(BUILDDIR)/Makefile $(BUILDDIR)/makefile_options $(LIBFILES) $(SRCDIR)/Makefile $(SRCDIR)/makefile_options
	make -C $(BUILDDIR) all VPATH=../src

gadget2_worker: worker_code.cc worker_code.h $(CODELIB) $(INTERFACE_OBJ) $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) $(SCINC) $(GSL_FLAGS) $< $(INTERFACE_OBJ) -o $@  $(CODELIB) $(SCLIBS) $(GSL_LIBS) -lgsl -lgslcblas $(ALLVARS)
	
gadget2_worker_sockets: worker_code-sockets.cc worker_code.h $(CODELIB) $(INTERFACE_OBJ) $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) $(SCINC) $(GSL_FLAGS) $< $(INTERFACE_OBJ) -o $@  $(CODELIB) $(SCLIBS) $(GSL_LIBS) -lgsl -lgslcblas $(ALLVARS)
	

$(PERIODIC_BUILDDIR)/Makefile: $(SRCDIR)/Makefile 
	-mkdir $(PERIODIC_BUILDDIR) 
	cp $(SRCDIR)/Makefile $(PERIODIC_BUILDDIR)/Makefile
	
$(PERIODIC_BUILDDIR)/makefile_options: $(SRCDIR)/makefile_options 
	-mkdir $(PERIODIC_BUILDDIR) 
	cp $(SRCDIR)/makefile_options $(PERIODIC_BUILDDIR)/makefile_options
	
$(PERIODIC_BUILDDIR)/libgadget.a:  $(PERIODIC_BUILDDIR)/Makefile $(PERIODIC_BUILDDIR)/makefile_options
	make -C $(PERIODIC_BUILDDIR) all OPT="$(PERIODIC_OPT)" VPATH=../src
	
gadget2_worker_periodic: worker_code.cc worker_code.h gadget_code.po $(PERIODIC_BUILDDIR)/libgadget.a  $(INTERFACE_SRC) $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) -DPERIODIC $(SCINC) $(GSL_FLAGS) $< gadget_code.po -o $@  $(PERIODIC_BUILDDIR)/libgadget.a  $(SCLIBS) $(GSL_LIBS) -lgsl -lgslcblas $(PERIODIC_BUILDDIR)/allvars.o

gadget2_worker_periodic_sockets: worker_code-sockets.cc worker_code.h gadget_code.po $(PERIODIC_BUILDDIR)/libgadget.a  $(INTERFACE_SRC) $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) -DPERIODIC $(SCINC) $(GSL_FLAGS) $< gadget_code.po -o $@  $(PERIODIC_BUILDDIR)/libgadget.a $(SCLIBS) $(GSL_LIBS) -lgsl -lgslcblas $(PERIODIC_BUILDDIR)/allvars.o
	
gadget_code.po: gadget_code.cc $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) $(SCINC) $(GSL_FLAGS) -DPERIODIC -DNOGRAVITY -c -o $@ $< $(SCLIBS) $(GSL_LIBS)
	
gadget_code.o: gadget_code.cc $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) $(SCINC) $(GSL_FLAGS) -c -o $@ $< $(SCLIBS) $(GSL_LIBS)

