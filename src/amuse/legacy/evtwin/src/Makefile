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
PYTHON_INCLUDE_DIR = /usr/include/python2.5

# Default FORTRAN compiler if FORT variable not defined. It is better to
# set FORT from the shell login script than to change the default here.
ifndef FORT
   FORT=ifort
endif

# Other tools
LD = ld
AR = ar ruv
RANLIB = ranlib

# Actually, just use the compiler to pass information to the loader. This
# works better because the compiler will also include the proper library
# dependancies.
LD = $(FORT)
LDFLAGS = -shared
LIBS = -lm

# C compiler options, for the SWIG C wrapper for the MUSE plugin
CC = gcc
CFLAGS = -O3 -fomit-frame-pointer -fPIC

# SWIG
SWIG = swig

# Determine appropriate compiler settings for different compilers
# Compiler flags for the PortlandGroup compiler
ifeq ($(notdir $(FORT)),pgf95)
   FFLAGS += -O2 -Mextend
   FFLAGS += -module $(MODDIR) -I$(MODDIR)
endif

# Compiler flags for the g95 compiler
ifeq ($(notdir $(FORT)),g95)
   FFLAGS += -ffixed-line-length-132 -I$(MODDIR) -fmod=$(MODDIR)
   FFLAGS += -fzero -O -finline -pipe
	FLIB = -lg95
endif

# Compiler flags for the gfortran compiler
ifeq ($(notdir $(FORT)),gfortran)
   FFLAGS += -ffixed-line-length-132 -I$(MODDIR) -J$(MODDIR)
   FFLAGS += -O2
	FLIB = -lgfortran

ifdef OPENMP
   FFLAGS += -fopenmp
endif
endif

# Compiler flags for the Lahey compiler
ifeq ($(notdir $(FORT)),lf95)
   FFLAGS += -O --wide -I$(MODDIR) -M$(MODDIR)
endif

# Compiler flags for the Intel compiler
ifeq ($(notdir $(FORT)),ifort)
   #FFLAGS += -fltconsistency
   FFLAGS += -fpconstant
   FFLAGS += -extend_source
   FFLAGS += -align all
   FFLAGS += -module $(MODDIR) -I$(MODDIR)
   FFLAGS += -ip

ifdef OPTIMISE_MORE
   FFLAGS += -O3
#   FFLAGS += -no-inline-factor
else
   FFLAGS += -O2
endif

   FFLAGS += -axSSE4.2,SSSE3,SSE2

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

ifdef DEBUGMODE
   FFLAGS += -g 
   ifeq ($(FORT),ifort)
      FFLAGS += -fpe0 -fpstkchk
   endif
endif

ifdef PROFILEMODE
   FFLAGS += -g -pg
endif

include makefile.sys

# Set architecture flag
FFLAGS += $(ARCHITECTURE_FLAG)
CFLAGS += $(ARCHITECTURE_FLAG)

# Main source
SHARED_SRC = initdat.f initrun.f nremesh.f extrapolate_dh.f mesh.f control.f xtraelem.f fcompare.f consts.f menc.f plotvars.f bbegin.f cfuncs.f estate.f massloss.f file_exists.f positron_fdirac.f printb.f fgb2hb_comp.f
SHARED_SRC += dsolve.f 
RAWSRC = aamain.f $(SHARED_SRC)
LIBSRC = library.f $(SHARED_SRC) prtoft.f import.f export.f

# Target object files
OBJS = $(addprefix $(OBJDIR)/,$(RAWSRC:.f=.o))
FPIC_OBJS = $(addprefix $(OBJDIR)/pic_,$(LIBSRC:.f=.o))
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

depend:
	@echo Generating dependencies...
	@perl depend.pl code/*.f > makefile.dep
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

include makefile.dep

.f.o:
	$(FORT) $(FFLAGS) -c $?

# Executable and library targets
$(EV): $(OBJS)
	$(FORT) $(FFLAGS) -o $@ $(OBJS)

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

distclean: clean

veryclean: distclean
	rm -f makefile.dep
	
