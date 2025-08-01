!> This module contains variables that can be used to control the amount of
!! fudging that goes on in the code. Different types of fudges can be allowed and
!! disallowed.
!<
module control
  use real_kind
  use mesh

  ! Enable or disable fudge control
  logical :: use_fudge_control = .true.

  ! Allow or disallow underrelaxation to improve convergence
  logical :: allow_underrelaxation = .true.

  ! Allow or disallow overrelaxation to improve convergence
  logical :: allow_overrelaxation = .true.

  ! Allow or disallow relaxation of the mass-loss/gain rate
  ! Only operates on CMI for now, so the default setting is disabled.
  logical :: allow_mdotrelaxation = .false.

  ! Allow or disallow interpolation of the mean molecular weight
  ! This only really affects thermohaline mixing
  logical :: allow_avmurelaxation = .false.

  ! Decide whether to use the mean molecular weight on the current, or
  ! previous timestep. Default is the previous timestep, which is more
  ! numerically stable.
  logical :: use_previous_mu = .true.

  ! Allow or disallow relaxation of the angular momentumadvection term.
  ! See below for explanation.
  logical :: allow_amadvrelaxation = .false.

  ! Allow the code to "recycle" the current timestep
  ! This means: reduce the timestep, then use those corrections as a
  ! basis for a larger timestep. This way the code may avoid reducing
  ! the timestep and gain a net speed-increase.
  logical :: recycle_timestep = .false.

  ! Use the line search algorithm for global convergence described in
  ! Numerical Recipes (2nd edition) section 9.7 to scale the
  ! corrections of the Newton iteration step. Slower and doesn't work
  ! in TWIN mode (yet).
  logical :: use_linesearch = .false.

  ! Fudge factor for entropy adjustment from the artificial energy term
  real(double) :: impose_entropy_factor = 1.0D-7

  ! Fudge factor for the artificial composition adjustment. This is
  !  gradually turned up to 1 by printb.
  ! 1.0 means apply full corrections at each timestep.
  real(double) :: impose_composition_factor = 1.0d-4

  ! Try fudging the luminosity by supressing (or not) the TdS/dt contribution
  ! 1.0 means this term is completely taken along, 0.0 means it is completely ignored
  real(double) :: luminosity_fudge = 1.0

  ! Try reducing the strength of convective mixing and then allow it to grow
  ! to its proper mixing-length model value.
  ! 1.0 is full strength.
  real(double) :: mixing_fudge = 1.0

  ! Supress (or not) problematic terms in the luminosity equation
  ! 1.0 uses the full luminosity equation, 0.0 reduces the problems
  real(double) :: llumi_smooth = 1.0

  ! Reduce the mass-loss/gain rate to make it easier to converge
  ! This only affects the CMI term at present
  ! 1.0 gives the full mass-loss/gain rate, 0.0 always keeps the mass constant
  real(double) :: mdot_smooth = 1.0

  ! Interpolate between the mean molecular weight at the current and
  ! previous timesteps. This smoothens the calculation of thermohaline
  ! mixing.
  ! Does a linear interpolation between the new value and the old value.
  ! 1.0 gives the new value, 0.0 gives the old value.
  real(double) :: avmu_smooth = 1.0

  ! Smoothly turn on the angular momentum advaction term
  real(double) :: adam_smooth = 1.0d0

  ! What convection scheme should we use?
  !   1: The original Eggleton 1972 scheme
  !   2: The Pols&Tout 2001 scheme
  integer :: convection_scheme = 1

  ! Should the molecular weight gradient be included in the convection
  ! criterion (Ledoux vs. Schwarzschild criterion)?
  real(double) :: convection_ledoux = 0.0d0

  ! What should CMI do?
  !   1: Make mass grow exponentially: mdot += m*CMI (ZAMS generation)
  !   2: Make mass grow linearly: mdot += CMI
  integer :: cmi_mode = 1

  ! Should we use a linear or a quadratic prediction for the changes in
  ! the independent variables between timesteps?
  logical :: use_quadratic_predictions = .true.

end module control
