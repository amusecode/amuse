!> This module contains variables that can be used to control the amount of 
!! fudging that goes on in the code. Different types of fudges can be allowed and
!! disallowed.
!<
module control
  use real_kind
  use mesh
  
  ! Enable or disable fudge control
  logical :: use_fudge_control = .true.
  
  ! Allow or disallow fudges to be used
  logical :: fudge_allowed = .true.
  
  ! Allow or disallow extension if close to convergence on the last iteration
  ! NO LONGER USED, VARIABLE RETAINED FOR BACKWARD COMPATIBILITY
  logical ::allow_extension = .true.
  
  ! Allow or disallow underrelaxation to improve convergence
  logical :: allow_underrelaxation = .true.
  
  ! Allow or disallow overrelaxation to improve convergence
  logical :: allow_overrelaxation = .true.
  
  ! Allow or disallow interpolation of energy generation
  ! NO LONGER USED, VARIABLE RETAINED FOR BACKWARD COMPATIBILITY
  logical :: allow_egenrelaxation = .true.
  
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
  
  ! Apply or don't apply second-order corrections to the moving mesh 
  ! advection term in the diffusion equations.
  logical :: apply_second_order_corrections = .false.
  
  ! Allow the code to "recycle" the current timestep
  ! This means: reduce the timestep, then use those corrections as a
  ! basis for a larger timestep. This way the code may avoid reducing
  ! the timestep and gain a net speed-increase.
  logical :: recycle_timestep = .false.
  
  ! Allow the thermal energy generation term to be constructed by
  ! interpolation from the previous model (at constant mass), which
  ! eliminates the advection term.
  logical :: interpolate_thermal_energy = .false.
  
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
  
  ! Try boosting the convective mixing to find the appropriate corrections
  ! for the next timestep. Does a linear interpolation between the proper
  ! minxing-length value and a boosted rate.
  ! 1.0 is the normal, non-boosted mixing length rate.
  real(double) :: mixing_boost = 1.0
  
  ! Try smoothening the transition of energy generation between timesteps
  ! Does a linear interpolation between the new value and the old value.
  ! 1.0 gives the new value, 0.0 gives the old value.
  real(double) :: egen_smooth = 1.0
  
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
  
  ! Smoothing term for advection term in angular momentum transport.
  ! Angular momentum transport becomes complicated (for instance) near
  ! the base of the giant branch: the core is contracting, but at the same
  ! time meshpoints are migrating from the core to the hydrogen burning
  ! shell. Both of these effects affect the advection term in the equation
  ! for angular momentum transport non-linearly. The problem is compounded
  ! if the region is rigidly rotating and there is no other mixing process
  ! operating, because then the first correction has dw/dk=0, but the second
  ! iteration has not, which activates the (dm/dt dw/dk) term on the second
  ! iteration. A "proper" fix would be to insert a better guess for dw/dt,
  ! but for now we will settle on relaxing the advection term itself.
  ! 1.0 gives the full value of the advection term, 0.0 ignores it.
  real(double) :: amadv_smooth = 1.0
  
  ! Sugimoto (1971) coefficient scaling pre factor. The ratio between the
  ! timestep and the local diffusion timescale is multiplied by this factor.
  ! Set it to something outrageously large (1e16 or thereabouts) to disable
  ! this feature. Don't try anything below 1.
  real(double) :: off_centre_weight = 1.0D16
  
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
  
  !     Use reformulation of boundary conditions for TWIN into more
  !     intuituve variable names. It is supposed to do the same as
  !     the old implementation. For spin up you will need the new boundary
  !     conditions. -SdM
  logical :: COMPUTE_TWINBC_THEOLDWAY=.TRUE.
  
end module control
