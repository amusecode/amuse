SHELL = /bin/sh
################################################################################
#
# User Options
################################################################################
FC =gfortran
#FC=ifort

f2pyFC=gnu95 #this must match FC (ifort=intelem, gfortran=gnu95)

#Unforgiving debugging flags
#FFLAGS =-g -fbacktrace -Wall -fcheck=all -fPIC
#Fast optimizing flags
FFLAGS = -O3 -fPIC -ffree-line-length-0


###############################################################################
EXT_SUFFIX = $(shell python3-config --extension-suffix)


PHYSICS_OBS= physics-core.o cloud.o hotcore.o sputtering.o cshock.o jshock.o collapse.o
CHEMISTRY_OBS = photoreactions.o surfacereactions.o chemistry.o
##simple makefile for uclchem
##user must point the FC variable to their preferred fortran compiler
##builds ode solver, physics module and chemistry module before linking together for main

##User can also compile a python module which contains the subroutines in wrap.f90 as functions
##to do this run "python -m numpy.f2py -c --help-fcompiler" to find your fortran compiler
##and edit the python makef

##physics module selected by changing physics variable to chosen fortran file.
main: constants.o  network.o ${PHYSICS_OBS} dvode.o ${CHEMISTRY_OBS} io.o wrap.f90 defaultparameters.f90
	${FC} ${FFLAGS} -o ../../uclchem io.o ${PHYSICS_OBS} dvode.o ${CHEMISTRY_OBS} network.o wrap.f90 main.f90

# The order here is the compilation order and it's vital that modules are compiled after their dependencies
python: constants.o  network.o ${PHYSICS_OBS} dvode.o ${CHEMISTRY_OBS} io.o wrap.f90  defaultparameters.f90
	python3 -m numpy.f2py -c --fcompiler=${f2pyFC} io.o ${PHYSICS_OBS} dvode.o ${CHEMISTRY_OBS} network.o -m uclchemwrap wrap.f90

chemistry.o: odes.f90 rates.f90 chemistry.f90
	${FC} ${FFLAGS} -c chemistry.f90

sputtering.o: surfacereactions.o sputtering.f90
	${FC} ${FFLAGS} -c sputtering.f90

%.o: %.f90 
	${FC} ${FFLAGS} -c $<

clean: 
	rm *.o *.mod *.so ../uclchem/*.so