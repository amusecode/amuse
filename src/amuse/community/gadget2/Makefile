# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include ${AMUSE_DIR}/config.mk

MPICXX ?= mpicxx

SRCDIR		= src
include $(SRCDIR)/makefile_options

BUILDDIR = build
PERIODIC_BUILDDIR = build_periodic
PERIODIC_NOGRAVITY_BUILDDIR = build_periodic_nogravity
NOGRAVITY_BUILDDIR = build_nogravity

PERIODIC_OPT := $(OPT) -DPERIODIC
PERIODIC_NOGRAVITY_OPT := $(OPT) -DPERIODIC -DNOGRAVITY
NOGRAVITY_OPT := $(OPT) -DNOGRAVITY

CXXFLAGS += -g -DTOOLBOX $(OPT)

CODE_GENERATOR = $(AMUSE_DIR)/build.py

#PROFLIBS ?= -L$(AMUSE_DIR)/lib/ibis/src/profiling -libisprofiling -lmpi

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

INTERFACE_OBJ = interface.o
INTERFACE_SRC  = $(INTERFACE_OBJ:.o=.cc)

all:	gadget2_worker


clean:
	rm -f *.o *.po *.pngo *.ngo *.pyc worker_code.cc worker_code.h worker_code worker_code_periodic
	rm -f gadget2_worker gadget2_worker_periodic gadget2_worker_sockets worker_code-sockets.cc
	rm -f gadget2_worker_periodic_sockets gadget2_worker_nogravity_sockets gadget2_worker_nogravity
	rm -f gadget2_worker_periodic_nogravity
	rm -rf *.dSYM
	make -C $(SRCDIR) clean
	rm -Rf $(PERIODIC_BUILDDIR)
	rm -Rf $(PERIODIC_NOGRAVITY_BUILDDIR)
	rm -Rf $(NOGRAVITY_BUILDDIR)
	rm -Rf $(BUILDDIR)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c $^ Gadget2Interface -o $@ 

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
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(GSL_FLAGS) $< $(INTERFACE_OBJ) -o $@  $(CODELIB) $(SC_MPI_CLIBS) $(GSL_LIBS) -lgsl -lgslcblas $(ALLVARS)
	
interface.o: interface.cc $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(GSL_FLAGS) -c -o $@ $< $(SC_MPI_CLIBS) $(GSL_LIBS)


$(PERIODIC_BUILDDIR)/Makefile: $(SRCDIR)/Makefile 
	-mkdir $(PERIODIC_BUILDDIR) 
	cp $(SRCDIR)/Makefile $(PERIODIC_BUILDDIR)/Makefile
	
$(PERIODIC_BUILDDIR)/makefile_options: $(SRCDIR)/makefile_options 
	-mkdir $(PERIODIC_BUILDDIR) 
	cp $(SRCDIR)/makefile_options $(PERIODIC_BUILDDIR)/makefile_options
	
$(PERIODIC_BUILDDIR)/libgadget.a:  $(PERIODIC_BUILDDIR)/Makefile $(PERIODIC_BUILDDIR)/makefile_options
	make -C $(PERIODIC_BUILDDIR) all OPT="$(PERIODIC_OPT)" VPATH=../src
	
gadget2_worker_periodic: worker_code.cc worker_code.h interface.po $(PERIODIC_BUILDDIR)/libgadget.a  $(INTERFACE_SRC) $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) -DPERIODIC $(SC_FLAGS) $(GSL_FLAGS) $< interface.po -o $@  $(PERIODIC_BUILDDIR)/libgadget.a  $(SC_MPI_CLIBS) $(GSL_LIBS) -lgsl -lgslcblas $(PERIODIC_BUILDDIR)/allvars.o

interface.po: interface.cc $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(GSL_FLAGS) -DPERIODIC -c -o $@ $< $(SC_MPI_CLIBS) $(GSL_LIBS)
	

$(PERIODIC_NOGRAVITY_BUILDDIR)/Makefile: $(SRCDIR)/Makefile 
	-mkdir $(PERIODIC_NOGRAVITY_BUILDDIR) 
	cp $(SRCDIR)/Makefile $(PERIODIC_NOGRAVITY_BUILDDIR)/Makefile
	
$(PERIODIC_NOGRAVITY_BUILDDIR)/makefile_options: $(SRCDIR)/makefile_options 
	-mkdir $(PERIODIC_NOGRAVITY_BUILDDIR) 
	cp $(SRCDIR)/makefile_options $(PERIODIC_NOGRAVITY_BUILDDIR)/makefile_options
	
$(PERIODIC_NOGRAVITY_BUILDDIR)/libgadget.a:  $(PERIODIC_NOGRAVITY_BUILDDIR)/Makefile $(PERIODIC_NOGRAVITY_BUILDDIR)/makefile_options
	make -C $(PERIODIC_NOGRAVITY_BUILDDIR) all OPT="$(PERIODIC_NOGRAVITY_OPT)" VPATH=../src
	
gadget2_worker_periodic_nogravity: worker_code.cc worker_code.h interface.pngo $(PERIODIC_NOGRAVITY_BUILDDIR)/libgadget.a  $(INTERFACE_SRC) $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) -DPERIODIC -DNOGRAVITY $(SC_FLAGS) $(GSL_FLAGS) $< interface.pngo -o $@  $(PERIODIC_NOGRAVITY_BUILDDIR)/libgadget.a  $(SC_MPI_CLIBS) $(GSL_LIBS) -lgsl -lgslcblas $(PERIODIC_NOGRAVITY_BUILDDIR)/allvars.o

interface.pngo: interface.cc $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(GSL_FLAGS) -DPERIODIC -DNOGRAVITY -c -o $@ $< $(SC_MPI_CLIBS) $(GSL_LIBS)
	

$(NOGRAVITY_BUILDDIR)/Makefile: $(SRCDIR)/Makefile 
	-mkdir $(NOGRAVITY_BUILDDIR)
	cp $(SRCDIR)/Makefile $(NOGRAVITY_BUILDDIR)/Makefile

$(NOGRAVITY_BUILDDIR)/makefile_options: $(SRCDIR)/makefile_options 
	-mkdir $(NOGRAVITY_BUILDDIR)
	cp $(SRCDIR)/makefile_options $(NOGRAVITY_BUILDDIR)/makefile_options

$(NOGRAVITY_BUILDDIR)/libgadget.a:  $(NOGRAVITY_BUILDDIR)/Makefile $(NOGRAVITY_BUILDDIR)/makefile_options
	make -C $(NOGRAVITY_BUILDDIR) all OPT="$(NOGRAVITY_OPT)" VPATH=../src

gadget2_worker_nogravity: worker_code.cc worker_code.h interface.ngo $(NOGRAVITY_BUILDDIR)/libgadget.a  $(INTERFACE_SRC) $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) -DNOGRAVITY $(SC_FLAGS) $(GSL_FLAGS) $< interface.ngo -o $@  $(NOGRAVITY_BUILDDIR)/libgadget.a  $(SC_MPI_CLIBS) $(GSL_LIBS) -lgsl -lgslcblas $(NOGRAVITY_BUILDDIR)/allvars.o

interface.ngo: interface.cc $(SRCDIR)/makefile_options
	$(MPICXX) $(CXXFLAGS) $(SC_FLAGS) $(GSL_FLAGS) -DNOGRAVITY -c -o $@ $< $(SC_MPI_CLIBS) $(GSL_LIBS)


