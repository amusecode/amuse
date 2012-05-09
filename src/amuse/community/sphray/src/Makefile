##############################################################################
##############################################################################
#                                                                            #
#                SPHRAY Makefile - Written by Gabriel Altay                  #
#                                                                            #
##############################################################################
##############################################################################


# Runtime Macros 
#=============================================================================
#
# These control the particle data structure and the code operation 
# in the case of the recombination macros.  The default is for all
# these flags to be commented out. 

#OPT += -DincHe      # if you want to include Helium
#OPT += -DoutGammaHI # if you want to output HI photoionization rate

#OPT += -DincHmf     # if you want individual Hydrogen mass fractions
#OPT += -DincHemf    # if you want individual Helium mass fractions

#OPT += -DincEOS     # store equation of state variable 
#OPT += -DincSFR     # store star formation rate
#OPT += -DincCloudy  # store CLOUDY eq. ionization values for particles

#OPT += -DuseHDF5    # if you are using the hdf5 library for anything

#OPT += -DOWLS       # selects input as an OWLS snapshot
#OPT += -DGIMIC      # selects input as a GIMIC snapshot


ifeq (OWLS,$(findstring OWLS, $(OPT)))	
   OPT += -DincHmf
   OPT += -DuseHDF5     
   OPT += -DoutGammaHI
   OPT += -DincCloudy
   OPT += -DincEOS
endif


ifeq (GIMIC,$(findstring GIMIC, $(OPT)))	
   OPT += -DincHmf
   OPT += -DuseHDF5     
   OPT += -DoutGammaHI
   OPT += -DincCloudy
   OPT += -DincEOS
endif



# System Definitions
#=============================================================================
#
# In order for the Makefile to work you have to design an include file
# for your system that sets certain variables.  Please look in the 
# makes directory for examples.  Right now, I have examples for three 
# fortran compilers.  gfortran is in make.blue-velvet.serial.opt, ifort
# is in make.santabuntu.serial.opt, and the sun compiler is in 
# make.icc-graphics.serial.opt
#

#include makes/make.icc-graphics.serial.debug
#include makes/make.icc-graphics.serial.opt

#include makes/make.xgpc.serial.debug
#include makes/make.xgpc.serial.opt

#include makes/make.ferrari.serial.opt

#include makes/make.cosma.serial.debug
#include makes/make.cosma.serial.opt

#include makes/make.santabuntu.serial.opt
#include makes/make.santabuntu.serial.debug

#include makes/make.blue-velvet.serial.gfortran.debug

#include makes/make.blue-velvet.serial.intel.debug
#include makes/make.blue-velvet.serial.intel.opt

#include makes/make.titania.serial.debug

#include makes/make.warp.opt

ifeq ($(findstring gfortran, $(notdir $(FORTRAN))), gfortran)
FC= gfortran
FFLAGS= -cpp -O0

ISLIB= -L
ISINC= -I
ISMOD= -I
ISRLIB= -L

DIRHDF= /home/galtay/usr/hdf5-1.8.5-patch1-gfortran-4.5/hdf5
endif


ifeq ($(findstring ifort, $(notdir $(FORTRAN))), ifort)
FC= ifort
IFORTFLAGS= -fpe0 -g -traceback -heap-arrays -fpp -u
FFLAGS= $(IFORTFLAGS) -O3 -no-prec-div


ISLIB= -L
ISINC= -I
ISMOD= -I
ISRLIB= -Wl,-rpath=

DIRHDF= /usr/lib64
endif


#=============================================================================
#=============================================================================
CC= gcc
CFLAGS= -O3


APPS= screen sphray 

#--------------------

SRC= 	myf03.o \
	sobol.o \
	physical_constants.o \
	mt19937.o \
	m_mrgrnk.o \
	hui_gnedin_atomic_rates.o \
	cen_atomic_rates.o \
	hummer_atomic_rates.o  \
	atomic_rates.o \
	ion_table_class.o \
        gadget_general_class.o \
	gadget_public_header_class.o \
        gadget_owls_header_class.o \
	gadget_public_header_hdf5_class.o \
        gadget_sphray_header_class.o \
	cosmology.o \
	particle_system.o \
	b2cd.o \
	spectra.o \
	sphpar.o \
	octtree3.o \
	ray.o \
	raylist.o \
	global.o \
	config.o \
	source_input.o \
	gadget_public_input.o \
        gadget_cosmoBH_input.o \
        gadget_owls_input.o \
        gadget_vbromm_input.o \
        gadget_public_input_hdf5.o \
        update_particles.o \
	main_input.o \
	ionpar.o \
	euler.o \
	bdf.o \
	iliev_comparison_project.o \
	output.o \
	ion_temperature_update.o \
	initialize.o \
	mainloop.o 


LIBFILES= $(SRC)  amuse_mainloop.o hash.o amuse_helpers.o

#=============================================================================
# These lines should NOT be edited if you are using the HDF5 libraries.
# Instead, simply set the variable DIRHDF, preferably in one of the 
# include files specified above.
#


ifeq (useHDF5,$(findstring useHDF5, $(OPT)))
   INCHDF= $(ISINC)$(DIRHDF)/include $(ISINC)$(DIRHDF)/lib
   MODHDF= $(ISMOD)$(DIRHDF)/include
   RLIBHDF= $(ISRLIB)$(DIRHDF)/lib
   LIBHDF= $(ISLIB)$(DIRHDF)/lib 
   OPTHDF= $(INCHDF) $(MODHDF) $(RLIBHDF) $(LIBHDF)  -lhdf5 -lhdfwrapper
endif
OPT += $(OPTHDF)

#=============================================================================
# These only need to be set if you are using the OpenGL visiualization tool
# Note that I have taken code from the f03gl project and placed it in a 
# directory called OpenGL within the SPHRAY distribution to save you the 
# trouble of downloading it. HOWEVER the work on that project has been done
# by 
#
# Anthony Stone <ajs1 at cam dot ac dot uk> and 
# Aleksandar Donev <aleks dot donev at gmail dot com>
#
#
# http://www-stone.ch.cam.ac.uk/pub/f03gl
#
#=============================================================================

# local OpenGL directory (should not need to be changed) 
#-----------------------------------------------------------------------------
LCLGLDIR = OpenGL

#  change this to correspond to your GLUT implementation
#-----------------------------------------------------------------------------
# GLUT= glut
# GLUT= openglut
  GLUT= freeglut

GLSRC=  $(LCLGLDIR)/GLUT_fonts.o \
	$(LCLGLDIR)/OpenGL_gl.o  \
	$(LCLGLDIR)/OpenGL_glu.o \
	$(LCLGLDIR)/OpenGL_${GLUT}.o 

THRDSRC= $(LCLGLDIR)/fthread.o $(LCLGLDIR)/pt.o $(LCLGLDIR)/ptf77.o 

VIEWSRC= $(GLSRC) $(THRDSRC)


#  Comment out for 32-bit architecture
#-----------------------------------------------------------------------------
#BITS= 64

#  The OpenGL / GLUT libraries 
#-----------------------------------------------------------------------------
DIROGL= ${ISLIB}/usr/lib${BITS} $(ISLIB)./$(LCLGLDIR) 
LIBOGL= -lglut -lGL -lGLU

#  The X11 libraries 
#-----------------------------------------------------------------------------
DIRX11= ${ISLIB}/usr/X11R6/lib${BITS}
LIBX11= -lXaw -lXt -lXmu -lXi -lXext -lX11 -lm -lXxf86vm


# multi threading and local GLUT library
#=============================================================================


DIRTHRD = ${ISLIB}./$(LCLGLDIR)
LIBTHRD = -lpthread -lpt

DIRSVIEW = $(DIROGL) $(DIRX11) $(DIRTHRD) 
LIBSVIEW = ${LIBOGL} ${LIBX11} $(LIBTHRD) -ldl


# Targets
#=============================================================================
#
# Implicit Rules
#----------------
#
# $@ = name of target 
# $< = name of first dependency
# $^ = name of all dependencies with duplicates removed
# $? = name of all dependencies newer than the target
# $+ = name of all dependencies w/o duplicates removed
#

F2OBJ= -c    # Fortran flag to compile to object without linking
FNAME= -o    # Fortran flag to name output file

C2OBJ= -c    # C flag to compile to object without linking
CNAME= -o    # C flag to name output file

all:$(APPS)


# Just Dummy Reporting
#=============================================================================
screen: Makefile
	@echo
	@echo "OPT=    " $(OPT)
	@echo "SYSTEM= " $(SYSTEM)
	@echo

#
# Main SPHRAY application
#=============================================================================
sphray: $(SRC) sphray.o
	$(FC) $(FFLAGS) $^ $(OPT) $(FNAME) $@



# OpenGL Targets
#=============================================================================


# Create libpt library (for starting a seperate thread)

$(LCLGLDIR)/libpt.a: $(THRDSRC) $(LCLGLDIR)/pt.h
	ar crv $(LCLGLDIR)/libpt.a $(LCLGLDIR)/pt.o $(LCLGLDIR)/ptf77.o
	ranlib $(LCLGLDIR)/libpt.a

# SPHRAY with OpenGL viewer

glsphray: $(SRC) $(VIEWSRC) viewer.o $(LCLGLDIR)/libpt.a glsphray.o  
	${FC} $(FFLAGS) $(OPT) $^ $(DIRSVIEW) $(LIBSVIEW) $(FNAME) $@


# HDF5 modules
#=============================================================================
gadget_input_hdf5.o: gadget_input_hdf5.F90
	$(FC) $(FFLAGS) $(OPT) $(F2OBJ) $< $(FNAME) $@

# Test density estimates (in progress) 
#=============================================================================

density_test: $(SRC) density_test.o
	$(FC) $(FFLAGS) $(OPT) -o $@ $^ 


#
# AMUSE

amuse_interface:  libsphray.a

libsphray.a: $(LIBFILES)
	ar crs $@ $^


# Implicit Rules
#=============================================================================

%.o: %.F90 
	$(FC) $(FFLAGS) $(OPT) $(F2OBJ) $< $(FNAME) $@

%.o: %.f90 
	$(FC) $(FFLAGS) $(F2OBJ) $< $(FNAME) $@

%.o: %.f 
	$(FC) $(FFLAGS) $(F2OBJ) $< $(FNAME) $@

%.o: %.c
	$(CC) $(CFLAGS) $(C2OBJ) $< $(CNAME) $@




# Standard Cleaning Targets
#=============================================================================

clean :
	rm -f *.o *.mod *.bck *__genmod.f90 \
	$(LCLGLDIR)/libpt.a $(LCLGLDIR)/*.o $(LCLGLDIR)/*.mod

cleanall :
	rm -f *.o *.mod *.bck *__genmod.f90 \
	$(LCLGLDIR)/libpt.a $(LCLGLDIR)/*.o $(LCLGLDIR)/*.mod \
	density_test glsphray libsphray.a $(APPS) 

tidy :
	rm -f *~ 












