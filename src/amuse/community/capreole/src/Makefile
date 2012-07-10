# Makefile for Capreole.
#
# Author: Garrelt Mellema

# This Makefile can make different versions of Capreole.
# These versions differ in their parallelization and/or
# different initial conditions (called "problems").
#
# Note 1: Parallelization
# The parallelization intended is specified in the name
# of the executable: _omp means OpenMP (shared memory), 
# _mpi means MPI (distributed memory). Both can also be
# used at the same time (if your architecture supports
# it.
#
# Note 2: Initial conditions
# Different initial conditions are specified by modules
# called "geometry"-"problem name", where geometry for
# 3D Capreole is almost always "cart" for cartesian.
# So, cart-densityfield is a problem module (initial
# conditions) for a triaxial cloud, possibly hit by
# a shock wave.
#
# Note 3: Compiler & Flags
# The compiler is specified by the FC variable (MPIFC for the MPI
# compiler). We have only extensively used the Intel F90 compiler. 
# Support for other compilers will have to be added.
# Parts of the code need to know about the compiler, this is
# done through preprocessor statements. So when compiling with
# intel compiler, -DIFORT needs to be specified. Support for
# new compilers thus needs to be added in the code too.
#
# Note 4: Recompiling
# Some dependencies are through module parameters, and thus
# not recognized by make. Best practise is to run "make clean"
# before running "make".
#-------------------------------------------------------

# Compiler: gfortran
#FC = gfortran # GNU compiler
#MPIFC = mpif90 # MPI compiler

# F90 options (gfortran)
#GFORTFLAGS = -O3 -DGFORT -DMPILOG
# Processor dependent optimization
#F90FLAGS1 = $(GFORTFLAGS) 

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
MPI_FLAGS = -DMPI # 
#MPI_FLAGS = -DMPI -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -openmp # For gfortran compiler

#-------------------------------------------------------
# Compiler: ifort (Intel) best tested
FC = ifort # Intel compiler

MPIFC ?= mpif90

# F90 options (ifort)
IFORTFLAGS = -O -g -traceback -ip -fpp1 -DIFORT
#IFORTFLAGS = -O3 -vec_report -u -fpe0 -ipo -DIFORT -shared-intel #-check all -traceback
#IFORTFLAGS = -O3 -vec_report -u -fpe0 -ipo -mcmodel=medium -shared-intel -DIFORT #-check all -traceback
# Processor dependent optimization
#F90FLAGS1 = $(IFORTFLAGS) 
F90FLAGS1 ?= -xP $(IFORTFLAGS) 
#F90FLAGS1 = -xO $(IFORTFLAGS) 
#F90FLAGS1 = -xT $(IFORTFLAGS) # Laptop 
#F90FLAGS1 = -xB $(IFORTFLAGS)

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
#MPI_FLAGS = -DMPI # 
#MPI_FLAGS = -DMPI -DMPILOG # Add more (MPI node) diagnostic output
OPENMP_FLAGS = -openmp # For Intel compiler

#-------------------------------------------------------

# Compiler: Sun
#(possible problems with constant definition. Cannot have sqrt in constant
# definition)
#FC = f95 # Sun compiler
#MPIFC = mpif90 # MPI compiler

# F90 options (ifort)
#SUNFLAGS = -O3 -DSUN
# Processor dependent optimization
#F90FLAGS1 = $(SUNFLAGS) 
#F90FLAGS1 = -xW $(SUNFLAGS) 

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
#MPI_FLAGS = -DMPI # 
#MPI_FLAGS = $(MPI_FLAGS) -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -openmp # For Sun compiler

#-------------------------------------------------------

# PGI compiler (not recently used/tested)
#FC = pf90
#MPIFC = mpif77
#MPIFC = mpif90

# F90 options (pgi)
#PGIFLAGS = -O3 -fast -DPGI
#F90FLAGS1 = -tp barcelona-64  $(PGIFLAGS) # ranger processors

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -DMPI 
#MPI_FLAGS = $(MPI_FLAGS) -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -mp 

#-------------------------------------------------------

# Absoft compilers (not recently used/tested)
#FC = f90
#ABSOFTF90FLAGS = -O3 -cpu:opteron #-DMPI
#ABSOFTF90FLAGS = -keep -Wv,-Pprocs,1 -O3 -cpu:opteron #-DMPI
LIBS = #-lU77

#-------------------------------------------------------

#LDR     = $(F90)

OPTIONS = $(F90FLAGS)

LDFLAGS = $(OPTIONS) #-L/afs/astro.su.se/pkg/intel/Compiler/11.1/056/lib/intel64/
LIBS = -lirc

GLINCLUDE = -I/usr/include/GL \
-I/data/pelupes/libraries/f90gl-1.2.11/include/GL/

GLLIB = -L/data/pelupes/libraries/f90gl-1.2.11/lib \
-lf90GLU -lf90GL -lf90glut  -lGLU -lGL -lglut 

X11LIB = -L/usr/X11R6/lib64 -lXaw -lXt -lXmu -lXi -lXext -lX11

THREADLIB = -L. -lpthread -lpt


#-------------------------------------------------------

# list of objects we're using

CONSTANTS = mathconstants.o cgsconstants.o cgsphotoconstants.o cgsastroconstants.o

COSMOPARMS = cosmoparms.o

PROT = protection2.o

LOF = lof-2.2.o

ROESOL = roesol-adv.o

INTEGRATE_TVD = integrate_tvd.o

TVDSOL = tvdsolver.o

INTEGRATE_VLFVS = integrate_vlfvs.o

# CART problems

CART-CONSTANT = cart-constant.o

CART-ELLIPSECLUMP = cart-ellipseclump.o

CART-DENSITYFIELD = cart-densityfield.o

CART-ENRIQUE = cart-enrique2.o

CART-HALO = cart-halo.o

CART-MINIHALO = cart-minihalo.o

# RT-3D
COOLING = cooling.o
HCOOLING = cooling_h.o
HMCOOLING = cooling_hm.o
RT-3D_BASIC = romberg.o tped.o radiation.o clumping.o doric.o thermal.o 
RT-3D_BASIC_ME = romberg.o tped.o radiation_me.o clumping.o doric.o thermal.o 
RT-3D = $(RT-3D_BASIC) sourceprops.o ionic.o
RT-3D_OMP = $(RT-3D_BASIC) sourceprops.o ionic_openmp.o
RT-3D-PP = $(RT-3D_BASIC) sourceprops_pp.o ionic_pp.o
RT-3D-PP_ME = $(RT-3D_BASIC_ME) sourceprops_pp.o ionic_pp.o
RT-3D-PP-NA = sourceprops_pp.o rad_evolve_planeparallel_noav.o $(RT-3D_BASIC)

# amuse interface

amuse_interface: F90=$(FC)
amuse_interface: F90FLAGS = $(F90FLAGS1)
amuse_interface: libcapreole.a

amuse_interface_mpi: F90=$(MPIFC)
amuse_interface_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
amuse_interface_mpi: libcapreole_mpi.a

libcapreole.a: precision.o file_admin.o string.o sizes.o noscaling.o $(CONSTANTS) abundances.o atomic.o \
 no_mpi.o clocks.o  amuse_mesh.o amuse_cart-coords.o \
 hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o \
 amuse_problem.o integrate-strang.o ah3out.o evolve.o amuse_helpers.o
	ar crs $@ $^

libcapreole_mpi.a: precision.o file_admin.o string.o sizes.o noscaling.o $(CONSTANTS) abundances.o atomic.o \
 amuse_mpi.o clocks.o  amuse_mesh.o amuse_cart-coords.o \
 hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o \
 amuse_problem.o integrate-strang.o ah3out.o evolve.o amuse_helpers.o
	ar crs $@ $^



# gl ----------------------------------------------------------------------


gl: F90=$(FC)
gl: F90FLAGS = $(F90FLAGS1)
gl: viewer.o libpt.a fthread.o

glcapreole: F90=$(FC)
glcapreole: F90FLAGS = $(F90FLAGS1)
glcapreole: gl precision.o file_admin.o string.o sizes.o noscaling.o $(CONSTANTS) abundances.o atomic.o no_mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-amuse.o integrate-strang.o  ah3out.o evolve.o viewer.o fthread.o glcapreole.o
	$(F90) $(OPTIONS) -o $@ file_admin.o string.o no_mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-amuse.o integrate-strang.o ah3out.o evolve.o viewer.o fthread.o glcapreole.o  $(LIBS) $(GLLIB) $(X11LIB) $(THREADLIB)


# amuse -------------------------------------------------------------------

cart-amuse: F90=$(FC)
cart-amuse: F90FLAGS = $(F90FLAGS1)
cart-amuse : precision.o file_admin.o string.o sizes.o noscaling.o $(CONSTANTS) abundances.o atomic.o no_mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-amuse.o integrate-strang.o  ah3out.o evolve.o capreole.o
	$(F90) $(OPTIONS) -o $@ file_admin.o string.o no_mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-amuse.o integrate-strang.o ah3out.o evolve.o capreole.o $(LIBS)

mpi_cart-amuse: F90=$(MPIFC)
mpi_cart-amuse: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
mpi_cart-amuse : precision.o file_admin.o string.o sizes.o noscaling.o $(CONSTANTS) abundances.o atomic.o mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-amuse.o integrate-strang.o  ah3out.o evolve.o capreole.o
	$(F90) $(OPTIONS) -o $@ file_admin.o string.o mpi.o clocks.o mesh.o cart-coords.o hydro.o time.o cart-routines.o $(PROT) boundary.o $(LOF) $(ROESOL) no_ionic.o cart-amuse.o integrate-strang.o ah3out.o evolve.o capreole.o $(LIBS)

clean:
	rm -f *.o *.mod *.l *.il *.vo libpt.a
	rm -f mpi_cart-amuse cart-amuse glcapreole 
	rm -f libcapreole.a libcapreole_mpi.a

.f.o:
	$(F90) -c $(OPTIONS) $<

.f90.o:
	$(F90) -c $(OPTIONS) $<

.F90.o:
	$(F90) -c $(OPTIONS) $<

f.mod:
	$(F90) -c $(OPTIONS) $<


viewer.o: viewer.f90
	$(FC) -fpp $(OPTIONS) -I./ $(GLINCLUDE) -c -o $@ $<

libpt.a: pt.c ptf77.c pt.h
	cc -c $(PIC) pt.c
	cc -c $(PIC) ptf77.c
	ar crv libpt.a pt.o ptf77.o
	ranlib libpt.a


.SUFFIXES: .f90 .F90 .mod .o




