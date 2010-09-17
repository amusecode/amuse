!> This module gets abused a lot because the common block that preceded
!! it (tvbles) was used for many different things.
!! It stores (among other things):
!!  * Global properties of the star and the binary
!!  * Information about *1 and *2, for TWIN boundary conditions
!!  * Binary information for the current timestep in non-TWIN mode
!!  * Information passed around between routines in printb
!!  * The timestep control parameter and the artificial energy term
!! These should all be grouped and moved to their own module/struct, as
!! appropriate.
!<


module test_variables
   use real_kind
   implicit none
   
   ! Stellar and binary properties
   real(double) :: dt         ! Current timestep [s]
   real(double) :: age        ! Current age [yr]
   real(double) :: mc(2)      ! Mass scale for the core (for mesh spacing function)
   real(double) :: hspn(2)    ! Spin angular momentum, star 1 and star 2
   real(double) :: rlf(2)     ! Roche lobe filling factor, ln R*/RL
   real(double) :: zet(2)     ! Wind mass loss from the system [1e33g/s]
   real(double) :: xit(2)     ! Mass transfer to companion, by RLOF and wind
   real(double) :: tn(2)      ! Nuclear timescale; set in printb [s]
   real(double) :: be(2)      ! Binding energy of the stellar envelope [?]
   real(double) :: bm         ! Total mass in the binary [1e33 g]
   real(double) :: bper       ! Binary period [days]

   ! Numerics
   real(double) :: enc        ! Artificial energy generation rate
   real(double) :: cdd        ! Timestep control parameter
   integer :: jhold     ! If jhold<3, then don't change the timestep

   ! Properties of the binary, non-TWIN mode
   real(double) :: t0         ! ??? age difference with primary, non-TWIN mode
   real(double) :: m0         ! mass of primary, non-TWIN mode
   real(double) :: mta        ! Time derivative of primary mass, non-TWIN mode
   real(double) :: om0        ! "Other" mass, non-TWIN mode
   real(double) :: omta       ! Time derivative of "Other" mass, non-TWIN mode
   real(double) :: a0         ! Separation, non-TWIN mode
   real(double) :: ata        ! Time derivative of orbital angular momentum, non-TWIN mode
   real(double) :: e0         ! Eccentricity at start of timestep, non-TWIN mode
   real(double) :: eta        ! Time derivative of eccentricity, non-TWIN mode
   real(double) :: om         ! "Other" mass, mass of companion (non-TWIN only?)

   ! Stuff for printb (doesn't belong here)
   real(double) :: bp         ! Poloidal component of magnetic field (?) [?]
   real(double) :: horb       ! Orbital angular momentum
   real(double) :: ro         ! ?
   real(double) :: tfr        !> \todo Timescale for tidal friction [s?] FIXME: for only one star?
   real(double) :: ra2        ! Alfven radius (squared) [1e11 cm]
   real(double) :: secc       ! Eccentricity (why???)
   real(double) :: wmh        ! Total amount of mass in H [Msun?]
   real(double) :: wmhe       ! Total amount of mass in He [Msun?]
   real(double) :: mh         ! Mass of H exhausted core [Msun?]
   real(double) :: mhe        ! Mass of He exhausted core [Msun?]
   real(double) :: mco        ! Mass of C/O exhausted core [Msun?]
   real(double) :: lh         ! Luminosity due to H burning
   real(double) :: lhe        ! Luminosity due to He burning
   real(double) :: lc         ! Luminosity due to C burning (and later)
   real(double) :: lnu        ! Neutrino luminosity
   real(double) :: lth        ! Thermal energy release
   real(double) :: mcb(8)     ! Mass boundaries of convection zones
   real(double) :: msb(6)     ! Mass boundaries of semiconvection zones
   real(double) :: rcb(8)     ! Radius boundaries of convection zones (for convective turnover time)
   real(double) :: tct(8)     ! Crossing time for each convection zone (for convective turnover time)
   
   !> \todo FIXME: Are the variables rs, vmg, sm, tc used?
   real(double) :: rs         ! FIXME: Unused?
   real(double) :: vmg        ! FIXME: Unused?
   real(double) :: sm         ! FIXME: Unused?
   real(double) :: tc(2)      ! FIXME: Unused?

   ! Backup of previous values, for when the timestep is reduced
   real(double) :: prev(81), pprev(81)
   integer :: jm2, jm1 

end module test_variables
