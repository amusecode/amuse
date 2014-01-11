MPICXX ?= mpicxx
MPICC ?= mpicc
GSL_FLAGS ?= $(shell gsl-config --cflags)
GSL_LIBS ?= $(shell gsl-config --libs)
GSL_INCL ?= $(GSL_FLAGS)

#----------------------------------------------------------------------
# From the list below, please activate/deactivate the options that     
# apply to your run. If you modify any of these options, make sure     
# that you recompile the whole code by typing "make clean; make".      
#                                                                      
# Look at end of file for a brief guide to the compile-time options.   
#----------------------------------------------------------------------

include makefile_options

# !!!                                                                                        !!!
# !!! Do not set the options below, but set them in the makefile_options file included above !!!
# !!!                                                                                        !!!

##--------------------------------------- Basic operation mode of code
##OPT   +=  -DPERIODIC
#OPT   +=  -DUNEQUALSOFTENINGS


##--------------------------------------- Things that are always recommended
#OPT   +=  -DPEANOHILBERT
#OPT   +=  -DWALLCLOCK


##--------------------------------------- TreePM Options
##OPT   +=  -DPMGRID=128
##OPT   +=  -DPLACEHIGHRESREGION=3
##OPT   +=  -DENLARGEREGION=1.2
##OPT   +=  -DASMTH=1.25
##OPT   +=  -DRCUT=4.5


##--------------------------------------- Single/Double Precision
##OPT   +=  -DDOUBLEPRECISION
##OPT   +=  -DDOUBLEPRECISION_FFTW


##--------------------------------------- Time integration options
#OPT   +=  -DSYNCHRONIZATION
##OPT   +=  -DFLEXSTEPS
##OPT   +=  -DPSEUDOSYMMETRIC
##OPT   +=  -DNOSTOP_WHEN_BELOW_MINTIMESTEP
##OPT   +=  -DNOPMSTEPADJUSTMENT


##--------------------------------------- Output 
##OPT   +=  -DHAVE_HDF5
##OPT   +=  -DOUTPUTPOTENTIAL
##OPT   +=  -DOUTPUTACCELERATION
##OPT   +=  -DOUTPUTCHANGEOFENTROPY
##OPT   +=  -DOUTPUTTIMESTEP


##--------------------------------------- Things for special behaviour
##OPT   +=  -DNOGRAVITY
##OPT   +=  -DNOTREERND
##OPT   +=  -DNOTYPEPREFIX_FFTW
##OPT   +=  -DLONG_X=60
##OPT   +=  -DLONG_Y=5
##OPT   +=  -DLONG_Z=0.2
##OPT   +=  -DTWODIMS
##OPT   +=  -DSPH_BND_PARTICLES
##OPT   +=  -DNOVISCOSITYLIMITER
##OPT   +=  -DCOMPUTE_POTENTIAL_ENERGY
##OPT   +=  -DLONGIDS
##OPT   +=  -DISOTHERM_EQS
##OPT   +=  -DADAPTIVE_GRAVSOFT_FORGAS
##OPT   +=  -DSELECTIVE_NO_GRAVITY=2+4+8+16

##--------------------------------------- Testing and Debugging options
##OPT   +=  -DFORCETEST=0.1


##--------------------------------------- Glass making
##OPT   +=  -DMAKEGLASS=262144


#----------------------------------------------------------------------
# Here, select compile environment for the target machine. This may need 
# adjustment, depending on your local system. Follow the examples to add
# additional target platforms, and to get things properly compiled.
#----------------------------------------------------------------------

#--------------------------------------- Select some defaults

CC       =  $(MPICC)  # sets the C-compiler
OPTIMIZE =  # sets optimization and warning flags
MPICHLIB =  # will be included by the linker -lmpich


#--------------------------------------- Select target computer

#SYSTYPE="MPA"
#SYSTYPE="Mako"
#SYSTYPE="Regatta"
#SYSTYPE="RZG_LinuxCluster"
#SYSTYPE="RZG_LinuxCluster-gcc"
#SYSTYPE="OpteronMPA"
#SYSTYPE="OPA-Cluster32"
#SYSTYPE="OPA-Cluster64"


#--------------------------------------- Adjust settings for target computer

ifeq ($(SYSTYPE),"MPA")
CC       =  mpicc   
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/usr/common/pdsoft/include
GSL_LIBS =  -L/usr/common/pdsoft/lib  -Wl,"-R /usr/common/pdsoft/lib"
FFTW_INCL= 
FFTW_LIBS= 
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  -lhdf5 -lz 
endif


ifeq ($(SYSTYPE),"OpteronMPA")
CC       =  mpicc   
OPTIMIZE =  -O3 -Wall -m64
GSL_INCL =  -L/usr/local/include
GSL_LIBS =  -L/usr/local/lib 
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL =  -I/opt/hdf5/include
HDF5LIB  =  -L/opt/hdf5/lib -lhdf5 -lz  -Wl,"-R /opt/hdf5/lib"
endif


ifeq ($(SYSTYPE),"OPA-Cluster32")
CC       =  mpicc   
OPTIMIZE =  -O3 -Wall 
GSL_INCL =  -I/afs/rzg/bc-b/vrs/opteron32/include
GSL_LIBS =  -L/afs/rzg/bc-b/vrs/opteron32/lib  -Wl,"-R /afs/rzg/bc-b/vrs/opteron32/lib"
FFTW_INCL=  -I/afs/rzg/bc-b/vrs/opteron32/include
FFTW_LIBS=  -L/afs/rzg/bc-b/vrs/opteron32/lib 
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  -lhdf5 -lz 
endif


ifeq ($(SYSTYPE),"OPA-Cluster64")
CC       =  mpicc   
OPTIMIZE =  -O3 -Wall -m64
GSL_INCL =  -I/afs/rzg/bc-b/vrs/opteron64/include
GSL_LIBS =  -L/afs/rzg/bc-b/vrs/opteron64/lib  -Wl,"-R /afs/rzg/bc-b/vrs/opteron64/lib"
FFTW_INCL=  -I/afs/rzg/bc-b/vrs/opteron64/include
FFTW_LIBS=  -L/afs/rzg/bc-b/vrs/opteron64/lib 
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  -lhdf5 -lz 
endif


ifeq ($(SYSTYPE),"Mako")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -O3 -march=athlon-mp  -mfpmath=sse
GSL_INCL = 
GSL_LIBS = 
FFTW_INCL= 
FFTW_LIBS= 
MPICHLIB =
HDF5INCL =  
HDF5LIB  =  -lhdf5 -lz 
endif


ifeq ($(SYSTYPE),"Regatta")
CC       =  mpcc_r 
OPTIMIZE =  -O5 -qstrict -qipa -q64
GSL_INCL =  -I/afs/rzg/u/vrs/gsl_psi64/include
GSL_LIBS =  -L/afs/rzg/u/vrs/gsl_psi64/lib                
FFTW_INCL=  -I/afs/rzg/u/vrs/fftw_psi64/include
FFTW_LIBS=  -L/afs/rzg/u/vrs/fftw_psi64/lib  -q64 -qipa
MPICHLIB =
HDF5INCL =  -I/afs/rzg/u/vrs/hdf5_psi64/include
HDF5LIB  =  -L/afs/rzg/u/vrs/hdf5_psi64/lib  -lhdf5 -lz 
endif


ifeq ($(SYSTYPE),"RZG_LinuxCluster")
CC       =  mpicci   
OPTIMIZE =  -O3 -ip     # Note: Don't use the "-rcd" optimization of Intel's compiler! (causes code crashes)
GSL_INCL =  -I/afs/rzg/u/vrs/gsl_linux/include
GSL_LIBS =  -L/afs/rzg/u/vrs/gsl_linux/lib  -Wl,"-R /afs/rzg/u/vrs/gsl_linux/lib"
FFTW_INCL=  -I/afs/rzg/u/vrs/fftw_linux/include
FFTW_LIBS=  -L/afs/rzg/u/vrs/fftw_linux/lib
HDF5INCL =  -I/afs/rzg/u/vrs/hdf5_linux/include
HDF5LIB  =  -L/afs/rzg/u/vrs/hdf5_linux/lib -lhdf5 -lz  -Wl,"-R /afs/rzg/u/vrs/hdf5_linux/lib"
endif


ifeq ($(SYSTYPE),"RZG_LinuxCluster-gcc")
CC       =  mpiccg
OPTIMIZE =  -Wall -g -O3 -march=pentium4
GSL_INCL =  -I/afs/rzg/u/vrs/gsl_linux_gcc3.2/include
GSL_LIBS =  -L/afs/rzg/u/vrs/gsl_linux_gcc3.2/lib  -Wl,"-R /afs/rzg/u/vrs/gsl_linux_gcc3.2/lib"
FFTW_INCL=  -I/afs/rzg/u/vrs/fftw_linux_gcc3.2/include
FFTW_LIBS=  -L/afs/rzg/u/vrs/fftw_linux_gcc3.2/lib  
HDF5INCL =  -I/afs/rzg/u/vrs/hdf5_linux/include
HDF5LIB  =  -L/afs/rzg/u/vrs/hdf5_linux/lib  -lhdf5 -lz  -Wl,"-R /afs/rzg/u/vrs/hdf5_linux/lib"
endif


ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(OPT))) 
HDF5INCL =
HDF5LIB  =
endif


OPTIONS =  $(OPTIMIZE) $(OPT)

EXEC   = Gadget2

LIBOBJS = run.o  predict.o begrun.o endrun.o global.o  \
	 init.o restart.o  io.o    \
	 accel.o   read_ic.o  ngb.o  \
	 system.o  allocate.o  density.o  \
	 gravtree.o hydra.o  driftfac.o  \
	 domain.o   potential.o  \
     forcetree.o   peano.o gravtree_forcetest.o \
	 pm_periodic.o pm_nonperiodic.o longrange.o \
	 density_at_pos.o \
	 timestep_limit.o  sigvel.o
#  ^^^ Durier & Dalla Vecchia timestepping fix ^^^

#	 timestep.o
#  ^^^ original timestepping ^^^

OBJS   = main.o allvars.o $(LIBOBJS)

INCL   = allvars.h  proto.h  tags.h  Makefile

CODELIB = libgadget.a

CFLAGS += $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL)


ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(OPT)))    # fftw installed with type prefix?
  FFTW_LIB = $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(OPT)))
  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif


#LIBS   =   $(HDF5LIB) -g  $(MPICHLIB)  $(GSL_LIBS) -lgsl -lgslcblas -lm $(FFTW_LIB) 
LIBS   =   $(HDF5LIB) -g  $(MPICHLIB)  $(GSL_LIBS) -lgsl -lgslcblas -lm

all:  $(EXEC) $(CODELIB)

$(EXEC): $(OBJS) makefile_options
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(CODELIB): $(LIBOBJS) makefile_options
	ar crs $@ $(LIBOBJS)

$(OBJS): $(INCL) makefile_options

clean:
	rm -f $(OBJS) $(EXEC) $(CODELIB)


#-----------------------------------------------------------------------
#
#   Brief guide to compile-time options of the code. More information
#   can be found in the code documentation.
#
# - PERIODIC:   
#     Set this if you want to have periodic boundary conditions.
#
# - UNEQUALSOFTENINGS:
#     Set this if you use particles with different gravitational
#     softening lengths.
#
# - PEANOHILBERT:    
#     This is a tuning option. When set, the code will bring the
#     particles after each domain decomposition into Peano-Hilbert
#     order. This improves cache utilization and performance.
#  
# - WALLCLOCK:       
#     If set, a wallclock timer is used by the code to measure internal
#     time consumption (see cpu-log file).  Otherwise, a timer that
#     measures consumed processor ticks is used.
#
# - PMGRID:     
#     This enables the TreePM method, i.e. the long-range force is
#     computed with a PM-algorithm, and the short range force with the
#     tree. The parameter has to be set to the size of the mesh that
#     should be used, (e.g. 64, 96, 128, etc). The mesh dimensions need
#     not necessarily be a power of two.  Note: If the simulation is
#     not in a periodic box, then a FFT method for vacuum boundaries is
#     employed, using an actual mesh with dimension twice(!) that
#     specified by PMGRID.
#
# - PLACEHIGHRESREGION: 
#     If this option is set (will only work together with PMGRID), then
#     the long range force is computed in two stages: One Fourier-grid
#     is used to cover the whole simulation volume, allowing the
#     computation of the longe-range force.  A second Fourier mesh is
#     placed on the region occupied by "high-resolution" particles,
#     allowing the computation of an intermediate scale force. Finally,
#     the force on short scales is computed with the tree. This
#     procedure can be useful for "zoom-simulations", provided the
#     majority of particles (the high-res particles) are occupying only
#     a small fraction of the volume. To activate this option, the
#     parameter needs to be set to an integer bit mask that encodes the
#     particle types that make up the high-res particles.
#     For example, if types 0, 1, and 4 form the high-res
#     particles, set the parameter to PLACEHIGHRESREGION=19, because
#     2^0 + 2^1 + 2^4 = 19. The spatial region covered by the high-res 
#     grid is determined automatically from the initial conditions. 
#     Note: If a periodic box is used, the high-res zone may not intersect
#     the box boundaries.
#
# - ENLARGEREGION: 
#     The spatial region covered by the high-res zone has a fixed size
#     during the simulation, which initially is set to the smallest
#     region that encompasses all high-res particles. Normally, the
#     simulation will be interrupted if high-res particles leave this
#     region in the course of the run. However, by setting this
#     parameter to a value larger than one, the size of the high-res
#     region can be expanded, providing a buffer region.  For example, 
#     setting it to 1.4 will enlarge its side-length by 40% (it remains
#     centered on the high-res particles). Hence, with this setting, the 
#     high-res region may expand or move by a limited amount. 
#     Note: If SYNCHRONIZATION is activated, the code will be able to
#     continue even if high-res particles leave the initial high-res
#     grid. In this case, the code will update the size and position of
#     the grid that is placed onto the high-resolution region
#     automatically. To prevent that this potentially happens every
#     single PM step, one should nevertheless assign a value slightly
#     larger than 1 to ENLARGEREGION.
#
# - ASMTH: 
#     This can be used to override the value assumed for the scale that
#     defines the long-range/short-range force-split in the TreePM
#     algorithm. The default value is 1.25, in mesh-cells.
#
# - RCUT: 
#     This can be used to override the maximum radius in which the
#     short-range tree-force is evaluated (in case the TreePM algorithm
#     is used). The default value is 4.5, given in mesh-cells.
#
# - DOUBLEPRECISION: 
#     This makes the code store and compute internal particle data in
#     double precision. Note that output files are nevertheless written
#     by converting the particle data to single precision.
#
# - DDOUBLEPRECISION_FFTW:
#     If this is set, the code will use the double-precision version of
#     FTTW, provided the latter has been explicitly installed with a
#     "d" prefix, and NOTYPEPREFIX_FFTW is not set. Otherwise the
#     single precision version ("s" prefix) is used.
#
# - SYNCHRONIZATION: 
#     When this is set, particles are kept in a binary hierarchy of
#     timesteps and may only increase their timestep if the new
#     timestep will put them into synchronization with the higher time
#     level.
#
# - FLEXSTEPS: 
#     This is an alternative to SYNCHRONIZATION. Particle timesteps are
#     here allowed to be integer multiples of the minimum timestep that
#     occurs among the particles, which in turn is rounded down to the
#     nearest power-of-two devision of the total simulated
#     timespan. This option distributes particles more evenly over
#     individual system timesteps, particularly once a simulation has
#     run for a while, and may then result in a reduction of work-load
#     imbalance losses.
#
# - PSEUDOSYMMETRIC: 
#     When this option is set, the code will try to "anticipate"
#     timestep changes by extrapolating the change of the acceleration
#     into the future. This can in certain idealized cases improve the
#     long-term integration behaviour of periodic orbits, but should
#     make little or no difference in most real-world applications. May
#     only be used together with SYNCHRONIZATION.
#
# - NOSTOP_WHEN_BELOW_MINTIMESTEP: 
#     If this is activated, the code will not terminate when the
#     timestep falls below the value of MinSizeTimestep specified in
#     the parameterfile. This is useful for runs where one wants to
#     enforce a constant timestep for all particles. This can be done
#     by activating this option, and by setting MinSizeTimestep and
#     MaxSizeTimestep to an equal value.
#
# - NOPMSTEPADJUSTMENT: 
#     When this is set, the long-range timestep for the PM-force
#     computation (when the TreePM algorithm is used) is always
#     determined by MaxSizeTimeStep.  Otherwise, it is determined by
#     the MaxRMSDisplacement parameter, or MaxSizeTimeStep, whichever
#     gives the smaller step.
#
# - HAVE_HDF5:
#     If this is set, the code will be compiled with support for input
#     and output in the HDF5 format. You need to have the HDF5
#     libraries and headers installed on your computer for this option
#     to work. The HDF5 format can then be selected as format "3" in
#     Gadget's parameterfile.
#
# - OUTPUTPOTENTIAL: 
#     This will make the code compute gravitational potentials for
#     all particles each time a snapshot file is generated. The values
#     are then included in the snapshot file. Note that the computation
#     of the values of the gravitational potential costs additional CPU.
#
# - OUTPUTACCELERATION: 
#     This will include the physical acceleration of each particle in
#     snapshot files.
#
# - OUTPUTCHANGEOFENTROPY: 
#     This will include the rate of change of entropy of gas particles
#     in snapshot files.
#
# - OUTPUTTIMESTEP:  
#     This will include the current timesteps of all particles in the 
#     snapshot files.
#
# - NOGRAVITY      
#     This switches off gravity. Useful only for pure SPH simulations 
#     in non-expanding space.
#
# - NOTREERND:       
#     If this is not set, the tree construction will succeed even when
#     there are a few particles at identical locations. This is done by
#     `rerouting' particles once the node-size has fallen below 1.0e-3
#     of the softening length. When this option is activated, this will
#     be surpressed and the tree construction will always fail if there
#     are particles at extremely close coordinates.
#
# - NOTYPEPREFIX_FFTW:
#     This is an option that signals that FFTW has been compiled
#     without the type-prefix option, i.e. no leading "d" or "s"
#     characters are used to access the library.
#
# - LONG_X/Y/Z: 
#     These options can be used together with PERIODIC and NOGRAVITY only.
#     When set, the options define numerical factors that can be used to
#     distorts the periodic simulation cube into a parallelepiped of 
#     arbitrary aspect ratio. This can be useful for idealized SPH tests.
#
# - TWODIMS:
#     This effectively switches of one dimension in SPH, i.e. the code
#     follows only 2d hydrodynamics in the xy-, yz-, or xz-plane. This
#     only works with NOGRAVITY, and if all coordinates of the third
#     axis are exactly equal. Can be useful for idealized SPH tests.
#
# - SPH_BND_PARTICLES:
#     If this is set, particles with a particle-ID equal to zero do not
#     receive any SPH acceleration. This can be useful for idealized
#     SPH tests, where these particles represent fixed "walls".
#
# - NOVISCOSITYLIMITER:   
#     If this is set, the code will not try to put an upper limit on
#     the viscous force in case an implausibly high pair-wise viscous
#     force (which may lead to a particle 'reflection' in case of poor
#     timestepping) should arise. Note: For proper settings of the
#     timestep parameters, this situation should not arise.
#
# - COMPUTE_POTENTIAL_ENERGY:
#     When this option is set, the code will compute the gravitational
#     potential energy each time a global statistics is computed. This
#     can be useful for testing global energy conservation.
#
# - LONGIDS: 
#     If this is set, the code assumes that particle-IDs are stored as
#     64-bit long integers. This is only really needed if you want to
#     go beyond ~2 billion particles.
#
# - ISOTHERM_EQS:
#     This special option makes the gas behave like an isothermal gas
#     with equation of state P = cs^2 * rho. The sound-speed cs is set by 
#     the thermal energy per unit mass in the intial conditions, 
#     i.e. cs^2=u. If the value for u is zero, then the initial gas 
#     temperature in the parameter file is used to define the sound speed
#     according to cs^2 = 3/2 kT/mp, where mp is the proton mass.
#
# - ADAPTIVE_GRAVSOFT_FORGAS:
#     When this option is set, the gravitational softening lengths used for
#     gas particles is tied to their SPH smoothing length. This can be useful
#     for dissipative collapse simulations. The option requires the setting
#     of UNEQUALSOFTENINGS.
#
# - SELECTIVE_NO_GRAVITY:
#     This can be used for special computations where one wants to 
#     exclude certain particle types from receiving gravitational
#     forces. The particle types that are excluded in this fashion are 
#     specified by a bit mask, in the same as for the PLACEHIGHRESREGION 
#     option.
#
# - FORCETEST:       
#     This can be set to check the force accuracy of the code. The
#     option needs to be set to a number between 0 and 1 (e.g. 0.01),
#     which is taken to specify a random fraction of particles for
#     which at each timestep forces by direct summation are
#     computed. The normal tree-forces and the correct direct
#     summation forces are collected in a file. Note that the
#     simulation itself is unaffected by this option, but it will of
#     course run much(!) slower, especially if 
#     FORCETEST*NumPart*NumPart >> NumPart. Note: Particle IDs must 
#     be set to numbers >=1 for this to work.
#
# - MAKEGLASS
#     This option can be used to generate a glass-like particle 
#     configuration. The value assigned gives the particle load,
#     which is initially generated as a Poisson sample and then
#     evolved towards a glass with the sign of gravity reversed.
#
#-----------------------------------------------------------------------

