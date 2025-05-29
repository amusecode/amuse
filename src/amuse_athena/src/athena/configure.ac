#-------------------------------------------------------------------------------
# configure.ac: Athena configuration script template file.  Run 'autoconf'
# to generate a new 'configure' script from this file.
#
# When configure is run, it uses the command line options and default settings
# to create custom versions of the following files:
#    Makefile, Makeoption, src/Makefile, src/defs.h, src/*/Makefile,
#
# If new configure options are added to this file, they should also be
# added in show_config() in src/show_config.c, and to the diagnostic message
# output at the end of this file.
#
# This template adds the following options to the 'configure' script:
#
# PHYSICS "packages":
#   --with-problem='name'                 (use prob/name.c as problem generator)
#   --with-gas=[hydro,mhd]                              (hydro or mhd algorithm)
#   --with-eos=[isothermal, adiabatic]                       (equation of state)
#   --with-nscalars=n                               (number of advected scalars)
#   --with-gravity=[fft,fft_obc,multigrid]                        (self-gravity)
#   --with-particles=[feedback,passive]              (dust particle integration)
#   --with-coord=[cartesian,cylindrical]                     (coordinate system)
#
# PHYSICS "features":
#   --enable-conduction                            (explicit thermal conduction)
#   --enable-resistivity                                  (explicit resistivity)
#   --enable-special-relativity              (special relativistic hydro or MHD)
#   --enable-viscosity                                      (explicit viscosity)
#
# ALGORITHM "packages":
#   --with-order=[1,2,3,2p,3p]        (order and type of spatial reconstruction)
#   --with-flux=[roe,hlle,hllc,hlld,force,exact,two-shock]       (flux function)
#   --with-integrator=[ctu,vl]                   (unsplit integration algorithm)
#   --with-cflags=[opt,debug,profile]                       (set compiler flags)
#
# ALGORITHM "features":
#   --enable-fargo                                      (enable FARGO algorithm)
#   --enable-fft                (compile and link with FFTW block decomposition)
#   --enable-fofc                 (first-order flux correction in VL integrator)
#   --enable-ghost                      (write out ghost cells in outputs/dumps)
#   --enable-h-correction              (turn on H-correction in multidimensions)
#   --enable-hllallwave                    (all-wave integration for HLL fluxes)
#   --enable-mpi                                          (parallelize with MPI)
#   --enable-shearing box                    (include shearing box source terms)
#   --enable-single                                 (double or single precision)
#   --enable-smr                                        (static mesh refinement)
#
#-------------------------------------------------------------------------------
# generic things

AC_INIT(src/athena.h)
AC_CONFIG_HEADER(src/config.h)
AC_PROG_CC
AC_PROG_INSTALL

#-------------------------------------------------------------------------------
# PHYSICS PACKAGE: name of problem generator from src/prob/name.c 
#   --with-problem=name (default is linear_wave1d)

AC_SUBST(PROBLEM)
AC_ARG_WITH(problem,
	[--with-problem=PROB  Problem generator from src/prob/name.c],
	with_problem=$withval, with_problem="linear_wave1d")
PROBLEM=$with_problem
if test -e src/prob/$PROBLEM.c; then
   rm -f src/problem.c
   ln -s prob/$PROBLEM.c src/problem.c
else
  ls src/prob
  AC_MSG_ERROR([Invalid problem filename, valid names listed above])
fi

#-------------------------------------------------------------------------------
# PHYSICS PACKAGE: gas mode (hydro or mhd)
#  --with-gas=hydro (default is mhd)
#   mode also sets the default flux handling, since that depends on whether the
#   problem is hydro or mhd (see option --with-flux= below)

AC_SUBST(GAS)
AC_ARG_WITH(gas,
	[--with-gas=GAS  Gas properties (mhd,hydro)],
	with_gas=$withval, with_gas="mhd")
if test   "$with_gas" = "hydro"; then
  GAS="HYDRO"
  flux_default="hllc"
elif test "$with_gas" = "mhd";   then
  GAS="MHD"
  flux_default="roe"
else
  AC_MSG_ERROR([expected --with-gas=mhd or hydro])
fi

#-------------------------------------------------------------------------------
# PHYSICS PACKAGE: equation of state
#   --with-eos=[adiabatic,isothermal] (default is adiabatic)

AC_SUBST(EOS)
AC_ARG_WITH(eos,
	[--with-eos=EOS  Equation of state (adiabatic,isothermal)],
	with_eos=$withval, with_eos=adiabatic)
if test   "$with_eos" = "adiabatic";  then
  EOS="ADIABATIC"
elif test "$with_eos" = "isothermal"; then
  EOS="ISOTHERMAL"
else
  AC_MSG_ERROR([expected --with-eos=adiabatic or isothermal])
fi

#-------------------------------------------------------------------------------
# PHYSICS PACKAGE: nscalars
#  --with-nscalars=n (n is any integer, default is 0)
#   number of passively advected scalars

AC_SUBST(NSCALARS)
AC_ARG_WITH(nscalars,
	[--with-nscalars=n (default is 0)],
	nscalars=$withval, nscalars=0)
if test   "$nscalars" = "0"; then
  NSCALARS="0"
elif test "$nscalars" != "0"; then
  NSCALARS=$nscalars
else
  AC_MSG_ERROR([expected --with-nscalars=n])
fi

#-------------------------------------------------------------------------------
# PHYSICS PACKAGE: self-gravity
#  --with-gravity= [fft,fft_obc,multigrid]

AC_SUBST(SELF_GRAVITY_DEFINE)
AC_SUBST(SELF_GRAVITY_ALGORITHM)
AC_ARG_WITH(gravity,
	[--with-gravity=SELF_GRAVITY_ALGORITHM  Algorithm for self gravity (fft, fft_obc, multigrid)],
	gravity_algorithm=$withval, gravity_algorithm=none)

if test   "$gravity_algorithm" = "fft"; then
  SELF_GRAVITY_DEFINE="SELF_GRAVITY"
  SELF_GRAVITY_USER="ON"
  SELF_GRAVITY_ALGORITHM="SELF_GRAVITY_USING_FFT"
elif test   "$gravity_algorithm" = "fft_obc"; then
  SELF_GRAVITY_DEFINE="SELF_GRAVITY"
  SELF_GRAVITY_USER="ON"
  SELF_GRAVITY_ALGORITHM="SELF_GRAVITY_USING_FFT_OBC"
elif test "$gravity_algorithm" = "multigrid"; then
  SELF_GRAVITY_DEFINE="SELF_GRAVITY"
  SELF_GRAVITY_USER="ON"
  SELF_GRAVITY_ALGORITHM="SELF_GRAVITY_USING_MULTIGRID"
elif test "$gravity_algorithm" = "none"; then
  SELF_GRAVITY_DEFINE="NO_SELF_GRAVITY"
  SELF_GRAVITY_USER="OFF"
  SELF_GRAVITY_ALGORITHM="SELF_GRAVITY_NONE"
else
  AC_MSG_ERROR([expected --with-gravity=fft, fft_obc, or multigrid])
fi

#-------------------------------------------------------------------------------
# PHYSICS PACKAGE: particles
#   --with-particles=[feedback,passive]

AC_SUBST(PARTICLES_ALGORITHM)
AC_SUBST(FEEDBACK_DEFINE)
AC_ARG_WITH(particles,
	[--with-particles=PARTICLES_ALGORITHM  particle integration algorithm (feedback, passive)], 
	particles_algorithm=$withval, particles_algorithm=none)

if test "$particles_algorithm" = "feedback"; then
  FEEDBACK_DEFINE="FEEDBACK"
  PARTICLES_USER="feedback"
  PARTICLES_ALGORITHM="PARTICLES"
elif test "$particles_algorithm" = "passive"; then
  FEEDBACK_DEFINE="NO_FEEDBACK"
  PARTICLES_USER="passive"
  PARTICLES_ALGORITHM="PARTICLES"
elif test "$particles_algorithm" = "none"; then
  FEEDBACK_DEFINE="NO_FEEDBACK"
  PARTICLES_USER="OFF"
  PARTICLES_ALGORITHM="NO_PARTICLES"
else
  AC_MSG_ERROR([expected --with-particles=feedback or passive])
fi

#-------------------------------------------------------------------------------
# PACKAGE: coordinate system
#   --with-coord=[cartesian,cylindrical] (default is cartesian)
  
AC_SUBST(COORD)
AC_ARG_WITH(coord,
	[--with-coord=COORD  Coordinate System (cartesian,cylindrical)],
	with_coord=$withval, with_coord=cartesian)
if test   "$with_coord" = "cartesian";  then
  COORD="CARTESIAN"
elif test "$with_coord" = "cylindrical"; then
  COORD="CYLINDRICAL"
else
  AC_MSG_ERROR([expected --with-coord=cartesian or cylindrical])
fi

#-------------------------------------------------------------------------------
# PHYSICS FEATURE: explicit thermal conduction
#  --enable-conduction
  
AC_SUBST(CONDUCTION_MODE) 
AC_ARG_ENABLE(conduction,
	[--enable-conduction  enable thermal conduction (default is no)],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  CONDUCTION_MODE="THERMAL_CONDUCTION"
  CONDUCTION_MODE_USER="ON"
else
  CONDUCTION_MODE="NO_THERMAL_CONDUCTION"
  CONDUCTION_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# PHYSICS FEATURE: explicit resistivity
#  --enable-resistivity
 
AC_SUBST(RESISTIVITY_MODE)
AC_ARG_ENABLE(resistivity,
	[--enable-resistivity  enable resistivity (default is no)],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  RESISTIVITY_MODE="RESISTIVITY"
  RESISTIVITY_MODE_USER="ON"
else
  RESISTIVITY_MODE="NO_RESISTIVITY"
  RESISTIVITY_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# PHYSICS FEATURE: special relativistic hydro or MHD
#  --enable-special-relativity

AC_SUBST(SPECIAL_RELATIVITY_MODE)
AC_ARG_ENABLE(special-relativity,
	[--enable-special-relativity  enable special relativistic hydro or MHD (default is no)],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  SPECIAL_RELATIVITY_MODE="SPECIAL_RELATIVITY"
  SPECIAL_RELATIVITY_MODE_USER="ON"
else
  SPECIAL_RELATIVITY_MODE="NEWTONIAN"
  SPECIAL_RELATIVITY_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# PHYSICS FEATURE: explicit viscosity
#  --enable-viscosity

AC_SUBST(VISCOSITY_MODE)
AC_ARG_ENABLE(viscosity,
	[--enable-viscosity  enable viscosity (default is no)],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  VISCOSITY_MODE="VISCOSITY"
  VISCOSITY_MODE_USER="ON"
else
  VISCOSITY_MODE="NO_VISCOSITY"
  VISCOSITY_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# ALGORITHM PACKAGE: set the order and type of the spatial reconstruction
#   --with-order=[1,2,3] - #-order reconstruction in the characteristic vars
#   --with-order=[2p,3p] - #-order reconstruction in the primitive      vars
#  (default is 2 -- 2nd-order reconstruction in the characteristic vars)

AC_SUBST(ORDER)
AC_SUBST(ACCURACY)

AC_ARG_WITH(order,
	[--with-order=ORDER      Order of Accuracy (1,2,3,2p,3p)],
	with_order=$withval, with_order="2")
ORDER=$with_order
if test   "$with_order" = "1"; then
  ACCURACY="FIRST_ORDER"
elif test "$with_order" = "2"; then
  ACCURACY="SECOND_ORDER_CHAR"
elif test "$with_order" = "2p"; then
  ACCURACY="SECOND_ORDER_PRIM"
elif test "$with_order" = "3"; then
  ACCURACY="THIRD_ORDER_CHAR"
elif test "$with_order" = "3p"; then
  ACCURACY="THIRD_ORDER_PRIM"
else
  AC_MSG_ERROR([expected --with-order=1, 2, 3, 2p, or 3p])
fi

#-------------------------------------------------------------------------------
# ALGORITHM PACKAGE: flux function
#   --with-flux=[roe,hlle,hllc,hlld,force,exact]
#   (default is roe if gas=mhd, hllc if gas=hydro)

AC_SUBST(FLUX_NAME)
AC_SUBST(FLUX_DEF)
AC_ARG_WITH(flux,
	[--with-flux=FLUX_NAME   Flux function (roe, hllc (hydro only), hlld, hlle, force, exact (hydro only),two-shock (hydro only))],
	with_flux=$withval, with_flux=$flux_default)
FLUX_NAME=$with_flux
if test   "$with_flux" = "roe"; then
  FLUX_DEF="ROE_FLUX"
elif test "$with_flux" = "hlle"; then
  FLUX_DEF="HLLE_FLUX"
elif test "$with_flux" = "hllc"; then
  FLUX_DEF="HLLC_FLUX"
elif test "$with_flux" = "hlld"; then
  FLUX_DEF="HLLD_FLUX"
elif test "$with_flux" = "force"; then
  FLUX_DEF="FORCE_FLUX"
elif test "$with_flux" = "exact"; then
  FLUX_DEF="EXACT_FLUX"
elif test "$with_flux" = "two-shock"; then
  FLUX_DEF="TWO_SHOCK_FLUX"
else
  AC_MSG_ERROR([Invalid flux function, valid types are: roe, hlle, hllc hlld, force, exact, two-shock])
fi

#-------------------------------------------------------------------------------
# ALGORITHM PACKAGE: set the unsplit integration algorithm
#   --with-integrator=[ctu,vl] (default is ctu)

AC_SUBST(INTEGRATOR)
AC_SUBST(INTEGRATOR_DEF)

AC_ARG_WITH(integrator,
	[--with-integrator=INTEGRATOR  Which unsplit integrator (ctu, vl)],
	with_integrator=$withval, with_integrator="ctu")
INTEGRATOR=$with_integrator
if test   "$with_integrator" = "ctu"; then
  INTEGRATOR_DEF="CTU_INTEGRATOR"
elif test "$with_integrator" = "vl"; then
  INTEGRATOR_DEF="VL_INTEGRATOR"
else
  AC_MSG_ERROR([expected --with-integrator=ctu or vl])
fi

#-------------------------------------------------------------------------------
# ALGORITHM PACKAGE: set compiler options.
#   --with-cflags=[opt,debug,profile] (default is opt)

AC_SUBST(DEBUG_MODE)
AC_SUBST(COMPILER_OPTS)

AC_ARG_WITH(cflags,
	[--with-cflags=DEBUG_MODE  Set debugging mode (opt,debug,profile)],
	with_debug=$withval, with_debug="opt")
if test "$with_debug" = "debug"; then
  COMPILER_OPTS="-g -Wall -W -ansi -pedantic"
  DEBUG_MODE="DEBUG"
elif test "$with_debug" = "profile"; then
  COMPILER_OPTS="-O3 -pg"
  DEBUG_MODE="PROFILE"
elif test "$with_debug" = "opt"; then
  COMPILER_OPTS="-O3"
  DEBUG_MODE="OPTIMIZE"
else
  AC_MSG_ERROR([expected --with-cflags=opt,debug, or profile])
fi


#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: precision of floating point arithmetic
#   --enable-single (default is double)

AC_SUBST(PRECISION)
AC_ARG_ENABLE(single,
	[--enable-single  single-precision (default is double)], 
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  PRECISION="SINGLE_PREC"
else
  PRECISION="DOUBLE_PREC"
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: write ghost cells in outputs/dumps
#   --enable-ghost

AC_SUBST(WRITE_GHOST_MODE)
AC_ARG_ENABLE(ghost,
	[--enable-ghost  write ghost zones],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  WRITE_GHOST_MODE="WRITE_GHOST_CELLS" 
  WRITE_GHOST_MODE_USER="ON" 
else
  WRITE_GHOST_MODE="NO_WRITE_GHOST_CELLS"
  WRITE_GHOST_MODE_USER="OFF" 
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: parallelize with MPI, --enable-mpi (default is no MPI)

AC_SUBST(MPI_MODE)
AC_ARG_ENABLE(mpi,
	[--enable-mpi  enable MPI parellelization],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  MPI_MODE="MPI_PARALLEL"
  MPI_MODE_USER="ON"
else
  MPI_MODE="NO_MPI_PARALLEL"
  MPI_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: turn on H-correction in multidimensional integrators
#   --enable-h-correction
# Note that H-correction only works with Roe flux.  Prints error message if
# H-correction is enabled with any other flux.

AC_SUBST(H_CORRECTION_MODE)
AC_ARG_ENABLE(h-correction,
	[--enable-h-correction  turn on H-correction],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
if test "$with_flux" = "roe"; then
  H_CORRECTION_MODE="H_CORRECTION"
  H_CORRECTION_MODE_USER="ON"
else
  AC_MSG_ERROR([H-correction only works with Roe flux])
fi
else
  H_CORRECTION_MODE="NO_H_CORRECTION"
  H_CORRECTION_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: turn on FFT libraries
#   --enable-fft

AC_SUBST(FFT_MODE)
AC_ARG_ENABLE(fft,
	[--enable-fft  compile and link FFT interface code (requires FFTW)],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  FFT_MODE="FFT_ENABLED"
  FFT_MODE_USER="ON"
else
  FFT_MODE="NO_FFT"
  FFT_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: turn on shearing box evolution
#   --enable-shearing-box
# Note shearing box only works with CTU+CT integrator.  Prints error message if
# shearing box is enabled with VL+CT integrator.

AC_SUBST(SHEARING_BOX_MODE)
AC_ARG_ENABLE(shearing-box,
	[--enable-shearing-box  turn on shearing-box],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
if test "$with_integrator" = "ctu"; then
  SHEARING_BOX_MODE="SHEARING_BOX"
  SHEARING_BOX_MODE_USER="ON"
else
  AC_MSG_ERROR([shearing-box only works with CTU+CT integrator])
fi
else
  SHEARING_BOX_MODE="NO_SHEARING_BOX"
  SHEARING_BOX_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: turn on FARGO algorithm
#   --enable-fargo

AC_SUBST(FARGO_MODE)
AC_ARG_ENABLE(fargo,
	[--enable-fargo  turn on fargo],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  FARGO_MODE="FARGO"
  FARGO_MODE_USER="ON"
else
  FARGO_MODE="NO_FARGO"
  FARGO_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: static mesh refinement
#   --enable-smr (default is no SMR)

AC_SUBST(MESH_REFINEMENT)
AC_ARG_ENABLE(smr,
	[--enable-smr  static mesh refinement (default is no refinement)], 
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  MESH_REFINEMENT="STATIC_MESH_REFINEMENT"
  SMR_MODE_USER="ON"
else
  MESH_REFINEMENT="NO_MESH_REFINEMENT"
  SMR_MODE_USER="OFF"
fi

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: turn on all-wave integration for HLL fluxes
#   --enable-hllallwaves
# This allows the L/R-states algorithms to integrate characteristics over
# domains of dependence of all waves, even those that lie on the "wrong" side of
# the interface.  Requires that one of the HLL-type lrstates algorithms be 
# enabled, otherwise prints error message.

AC_SUBST(HLL_ALL_WAVE_MODE)
AC_ARG_ENABLE(hllallwave,
	[--enable-hllallwave  turn on all-wave integration for HLL fluxes],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
if test "$with_flux" = "hlle"; then
  HLL_ALL_WAVE_MODE="HLL_ALL_WAVE"
  HLL_ALL_WAVE_MODE_USER="ON"
elif test "$with_flux" = "hllc"; then
  HLL_ALL_WAVE_MODE="HLL_ALL_WAVE"
  HLL_ALL_WAVE_MODE_USER="ON"
elif test "$with_flux" = "hlld"; then
  HLL_ALL_WAVE_MODE="HLL_ALL_WAVE"
  HLL_ALL_WAVE_MODE_USER="ON"
else
  AC_MSG_ERROR([All-wave integration only works with HLL-type fluxes!])
fi
else
  HLL_ALL_WAVE_MODE="NOT_HLL_ALL_WAVE"
  HLL_ALL_WAVE_MODE_USER="OFF"
fi  

#-------------------------------------------------------------------------------
# ALGORITHM FEATURE: first-order flux correction with VL integrator
#   --enable-fofc

AC_SUBST(FOFC_MODE)
AC_ARG_ENABLE(fofc,
	[--enable-fofc  enable first-order flux correction],
	ok=$enableval, ok=no)
if test "$ok" = "yes"; then
if test "$with_integrator" = "vl"; then
  FOFC_MODE="FIRST_ORDER_FLUX_CORRECTION"
  FOFC_MODE_USER="ON"
else
  AC_MSG_ERROR([first-order flux correction only works with VL integrator!])
fi
else
  FOFC_MODE="NO_FIRST_ORDER_FLUX_CORRECTION"
  FOFC_MODE_USER="OFF"
fi  

#-------------------------------------------------------------------------------
# check for compatibility of various options

if test "$with_coord" = "cylindrical"; then
  if test "$nscalars" != "0"; then
    AC_MSG_ERROR([Sorry, passive scalars and cylindrical coordinates are currently incompatible!])
  elif test "$gravity_algorithm" != "none"; then
    AC_MSG_ERROR([Sorry, self-gravity and cylindrical coordinates are currently incompatible!])
  elif test "$particles_algorithm" != "none"; then
    AC_MSG_ERROR([Sorry, particles and cylindrical coordinates are currently incompatible!])
  fi
fi

#-------------------------------------------------------------------------------
# check for various library functions

AC_PROG_GCC_TRADITIONAL
AC_CHECK_FUNCS(strdup)
AC_C_BIGENDIAN

#-------------------------------------------------------------------------------
# date and time of configure, becomes a macro in defs.h

A_CONFIGURE_DATE="`date`"
AC_SUBST(A_CONFIGURE_DATE)

#-------------------------------------------------------------------------------
# write final diagnostic output

AC_SUBST(WARNING1)
AC_SUBST(WARNING2)
WARNING1="WARNING! This file has been automatically generated by configure."
WARNING2="Any changes to it will be overwritten the next time configure is run."

AC_CONFIG_FILES(Makefile Makeoptions src/Makefile src/defs.h src/fftsrc/Makefile src/gravity/Makefile src/integrators/Makefile src/microphysics/Makefile src/reconstruction/Makefile src/rsolvers/Makefile src/particles/Makefile)
AC_OUTPUT

echo ""
echo "Your Athena distribution has now been configured:"
echo ""
echo "--- Physics Options ---"
echo ""
echo "Problem:                 $PROBLEM"
echo "Gas properties:          $GAS"
echo "Equation of State:       $EOS"
echo "Coordinate System:       $COORD"
echo "Advected scalar fields:  $NSCALARS"
if test "$gravity_algorithm" = "none"; then
  echo "Self-gravity:            $SELF_GRAVITY_USER"
else
  echo "Self-gravity:            $SELF_GRAVITY_ALGORITHM"
fi
echo "Resistivity:             $RESISTIVITY_MODE_USER"
echo "Viscosity:               $VISCOSITY_MODE_USER"
echo "Thermal conduction:      $CONDUCTION_MODE_USER"
echo "Particles:               $PARTICLES_USER"
echo "Special Relativity:      $SPECIAL_RELATIVITY_MODE_USER"
echo ""
echo "--- Algorithm Options ----"
echo ""
echo "Spatial Order:           $ORDER ($ACCURACY)"
echo "Flux:                    $FLUX_NAME"
echo "unsplit integrator:      $INTEGRATOR"
echo "Precision:               $PRECISION"
echo "Compiler options:        $COMPILER_OPTS"
echo "Ghost cell output:       $WRITE_GHOST_MODE_USER"
echo "Parallel modes: MPI      $MPI_MODE_USER"
echo "H-correction:            $H_CORRECTION_MODE_USER"
echo "FFT:                     $FFT_MODE_USER"
echo "Shearing-box:            $SHEARING_BOX_MODE_USER"
echo "FARGO:                   $FARGO_MODE_USER"
echo "All-wave integration:    $HLL_ALL_WAVE_MODE_USER"
echo "Static Mesh Refinement:  $SMR_MODE_USER"
echo "first-order flux corr:   $FOFC_MODE_USER"

