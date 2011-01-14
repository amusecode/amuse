.PRECIOUS: %.o

.PHONY: all libs dirs clean distclean veryclean compress depend old-depend install runtwin museplugin

# Directories for the object files, executable and generated module files
LIBDIR = lib
OBJDIR = obj
EXEDIR = code
MODDIR = modules

# Directory holding the source
SRCDIR = code

SONAME = so

# Destination directory for the MUSE plugin
TWINMUSEDIR = muse
PYTHON_INCLUDE_DIR = 

# Other tools
LD = ld
AR = ar ruv
RANLIB = ranlib

# C compiler options, for the SWIG C wrapper for the MUSE plugin
CC = gcc
CFLAGS = -O3 -fomit-frame-pointer -fPIC

# SWIG
SWIG = swig

include makefile.sys

# Default FORTRAN compiler if FORT variable not defined. It is best
# set FORT from the shell login script.
ifndef FORT
   FORT=$(FORTRAN)
endif

# Determine appropriate compiler settings for different compilers
# Compiler flags for the PortlandGroup compiler
ifeq ($(FORT),pgf95)
   FFLAGS += -O2 -Mextend
   FFLAGS += -module $(MODDIR) -I$(MODDIR)
endif

# Compiler flags for the g95 compiler
ifeq ($(FORT),g95)
   FFLAGS += -ffixed-line-length-132 -I$(MODDIR) -fmod=$(MODDIR)
   FFLAGS += -fzero -O -finline -pipe
	FLIB = -lg95
endif

# Compiler flags for the gfortran compiler
ifeq ($(FORT),gfortran)
   #FFLAGS += -ffixed-line-length-132 
   FFLAGS += -I$(MODDIR) -J$(MODDIR)
   FFLAGS += -O2 -Wall  
   FLIB = -lgfortran

ifdef OPENMP
   FFLAGS += -fopenmp
endif
endif

# Compiler flags for the Lahey compiler
ifeq ($(FORT),lf95)
   FFLAGS += -O --wide -I$(MODDIR) -M$(MODDIR)
endif


# Compiler flags for the Intel compiler
ifeq ($(FORT),ifort)
   #FFLAGS += -fltconsistency
   FFLAGS += -fpconstant
   FFLAGS += -extend_source
   FFLAGS += -align all
   FFLAGS += -module $(MODDIR) -I$(MODDIR)
   FFLAGS += -ip
   FFLAGS += -traceback

   ifdef OPTIMISE_MORE
      FFLAGS += -O3
      ifeq ($(FORT), ifort)
         FFLAGS += -unroll-agressive -vec-guard-write
      endif
#      FFLAGS += -no-inline-factor
   else

      ifdef RUNTIME_CHECKS
         FFLAGS += -O0 -warn all -check all -check noarg_temp_created
      else
         FFLAGS += -O2 -warn all  
      endif

   endif

   ifdef AS_FAST_AS_POSSIBLE
      #FFLAGS += -fp-model fast=2 -no-prec-div -fast-transcendentals #-fast
      FFLAGS += -no-prec-div -no-prec-sqrt -fast-transcendentals -fast
   else
      FFLAGS += -fp-model strict
   endif

   ifdef GLOBAL_OPTIMISATION
      FFLAGS += -ipo
   endif

   ifdef PARALLEL
      FFLAGS += -parallel
   endif

   ifdef STATICLINK
      FFLAGS += -static
   endif

   #FLIB = -lifcoremt -lifcore -lguide
   #FLIB = -limf -lifcore -lifport -lguide
endif
#End ifort


# Compiler flags for the NAG compiler
ifeq ($(FORT),nagfor)
   FFLAGS += -O0 -I$(MODDIR) -mdir $(MODDIR)
endif



ifdef DEBUGMODE
   FFLAGS += -g 
   ifeq ($(FORT),ifort)
      FFLAGS += -fpe0 -fpstkchk
   endif
   ifeq ($(FORT),gfortran)
      FFLAGS += -ffpe-trap=zero,invalid -fsignaling-nans
   endif
endif

ifdef PROFILEMODE
   FFLAGS += -g -pg
endif

# Set architecture flag
FFLAGS += $(ARCHITECTURE_FLAG)
CFLAGS += $(ARCHITECTURE_FLAG)

# Main source
SHARED_SRC = initdat.f90 initrun.f90 nremesh.f90 extrapolate_dh.f90 mesh.f90 control.f90 extra_elements.f90 fcompare.f90 constants.f90 
SHARED_SRC += menc.f90 ms2bss.f90 plotvariables.f90 begin.f90 update.f90 output_input_model.f90 stars_structure.f90 remesh.f90 eostate.f90 
SHARED_SRC += opacity.f90 opacity_co.f90 massloss.f90 file_exists.f90 positron_fdirac.f90 output_properties.f90 fgb2hb.f90 fgb2hb_comp.f90 
SHARED_SRC += nucleosynthesis.f90 semi_implicit.f90 interpolate.f90 diffusion.f90 step.f90 extra_solver.f90 eostate_cache.f90 pruner.f90 
SHARED_SRC += lt2ubv.f90 load_star_model.f90 svn_version.f90 neutrinos.f90 settings.f90 diffusion_coefficients.f90 test_variables.f90 
SHARED_SRC += atomic_data.f90 reaction_rate_data.f90 eostate_types.f90 funcs1_interface.f90 current_model_properties.f90 structure_variables.f90 
SHARED_SRC += binary_history.f90 accretion_abundances.f90 starting_values.f90 real_kind.f90 solver.f90
RAWSRC = ev_main.f90 $(SHARED_SRC)
LIBSRC = library.f90 $(SHARED_SRC) prtoft.f90 import.f90 export.f90

SRC = $(addprefix $(SRCDIR)/, $(LIBSRC) ev_main.f90)

# Target object files
OBJS1 = $(addprefix $(OBJDIR)/,$(RAWSRC:.f=.o ))
OBJS = $(OBJS1:.f90=.o)			# Hack
FPIC_OBJS1 = $(addprefix $(OBJDIR)/pic_,$(LIBSRC:.f=.o))
FPIC_OBJS = $(FPIC_OBJS1:.f90=.o)	# Hack
TWINMUSEOBJS = $(OBJDIR)/twinmuse.o $(OBJDIR)/twinmuse_wrap.o

EV = $(EXEDIR)/ev
LIBTWIN = $(LIBDIR)/libtwin.a
LIBTWINSO = $(LIBDIR)/libtwin.$(SONAME)
MUSEPLUGINLIB = _twinmuse.so

LIBS += $(FLIB) -lm -lpthread

all: makefile.dep dirs $(EV)

libs: $(LIBDIR) $(LIBTWIN) $(LIBTWINSO)

runtwin: $(OBJDIR)/twin.o libs
	$(FORT) -o $@ $< -Llib -ltwin $(FLIB)
#	$(FORT) -o $@ $< -static -Llib -ltwin

makefile.dep: Makefile
	@echo "Please run 'make depend' before running make"
	@exit

depend: $(SRCDIR)/svn_version.f90
	@echo Generating dependencies...
	@perl depend.pl $(SRC) > makefile.dep
	@echo Done!

DEPRULE = $(foreach src, $(RAWSRC), "$(addprefix $(OBJDIR)/,$(src:.f=.o)):$(addprefix $(SRCDIR)/,$(src))\n\t\$$(FORT) \$$(FFLAGS) -c $(addprefix \$$(SRCDIR)/,$(src)) -o $(addprefix \$$(OBJDIR)/,$(src:.f=.o))\n\n")
old-depend:
	@echo Generating dependencies...
	@echo -e $(DEPRULE) > makefile.dep
	@echo Done!

dirs: $(OBJDIR) $(LIBDIR) $(EXEDIR) $(MODDIR) $(TWINMUSEDIR)

$(TWINMUSEDIR):
	mkdir $@

$(LIBDIR):
	mkdir $@

$(OBJDIR):
	mkdir $@

$(EXEDIR):
	mkdir $@


$(MODDIR):
	mkdir $@

$(HOME)/bin/:
	mkdir $@

$(SRCDIR)/svn_version.f90:
	@perl detect_svn_version.pl

include makefile.dep

# Executable and library targets
$(EV): $(OBJS)
	$(FORT) $(FFLAGS) -o $@ $(OBJS)
	@rm $(SRCDIR)/svn_version.f90

$(LIBTWIN): $(FPIC_OBJS)
	$(AR) $@ $(FPIC_OBJS)

$(LIBTWINSO): $(FPIC_OBJS)
	$(LD) $(LDFLAGS) -o $@ $(FPIC_OBJS)

mkmergermod: $(OBJDIR)/mkmergermod.o $(LIBTWIN)
	$(FORT) $(FFLAGS) -L$(LIBDIR) -o $@ $(OBJDIR)/mkmergermod.o -ltwin $(LIBS)

install: all $(HOME)/bin/
	install $(EV) $(HOME)/bin/

# MUSE plugin
makefile.sys: configure.pl
	@echo "Please run ./configure.pl"
	@false

code/muse/config.h: $(LIBTWINSO) configure_twinmuse_wrapper
	./configure_twinmuse_wrapper
	mv config.h $@

code/muse/twinmuse_wrap.c:code/muse/config.h code/muse/twinmuse.c code/muse/twinmuse.i
	$(SWIG) -python code/muse/twinmuse.i

$(OBJDIR)/twinmuse_wrap.o:code/muse/twinmuse_wrap.c
	$(CC) $(CFLAGS) -c $< -o $@ -I$(PYTHON_INCLUDE_DIR)

$(OBJDIR)/twinmuse.o:code/muse/twinmuse.c code/muse/config.h
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBDIR)/$(MUSEPLUGINLIB): $(TWINMUSEOBJS)
	$(LD)  $(LDFLAGS) -o $@ $(FPIC_OBJS) $(TWINMUSEOBJS) $(LIBS)
#	$(CC) -shared -o $@ $(FPIC_OBJS) $(TWINMUSEOBJS) $(FLIB)

museplugin: $(TWINMUSEDIR) $(LIBTWINSO) $(LIBTWIN) $(LIBDIR)/$(MUSEPLUGINLIB)
	ln -sf $(PWD)/$(LIBDIR)/$(MUSEPLUGINLIB) $(TWINMUSEDIR)
	cp code/muse/twinmuse.py $(TWINMUSEDIR)

# Cleanup targets
clean:
	rm -f $(LIBDIR)/*.a
	rm -f $(LIBDIR)/*.so
	rm -f $(OBJDIR)/*.o
	rm -f $(TWINMUSEDIR)/*.pyc $(TWINMUSEDIR)/twinmuse.py $(TWINMUSEDIR)/*.so
	rm -f code/muse/config.h code/muse/twinmuse.py code/muse/twinmuse_wrap.c
	rm -f $(MODDIR)/*.mod
	rm -f core core.*
	rm -f $(EV)
	rm -f $(SRCDIR)/svn_version.f90

distclean: clean

veryclean: distclean
	rm -f makefile.dep

