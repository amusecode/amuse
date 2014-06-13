module stopping_conditions
use real_kind
implicit none

! Physical limits that trigger stopping conditions or special behaviour, from init.run
!
! uc(1)  - rlf1     Stop when *1 exceeds Roche lobe by this amount (measured as log(R*/RL))
! uc(2)  - age      Stop when age reaches this limit (in yr)
! uc(3)  - LCarb    Stop if luminosity due to C burning exceeds this limit and...
! uc(4)  - rlf2     ... the central density exceeds this limit (degenerate C ignition/C flash)
! uc(5)  - LHe      Start ZAHB construction if luminosity due to He burning exceeds this limit and...
! uc(6)  - rho      ... the central density exceeds this limit (degenerate He ignition/He flash)
! uc(7)  - MCO      Stop if the mass of the CO core exceeds this limit and...
! uc(8)  - rho      ... the central density exceeds this limit (degenerate core/core collapse?)
! uc(9)  - mdot     Stop if the mass loss rate exceeds this limit (units?)
! uc(10) - XHe      Unused
! uc(11) - He-eps   Unused 
! uc(12) - dtmin    Stop if the timestep becomes smaller than this limit (in s)
! uc(13) - sm8      Switch off artificial accretion (CMI) when exceeding this mass (for ZAHB construction)
! uc(14) - vmh8     Resume normal evolution (after ZAHB construction) once the core mass exceeds this mass
! uc(15) - XH       Stop the evolution when core H drops below this limit (if > 0)
! uc(16) - Rmax     Stop the evolution when the radius exceeds this limit (if > 0)
! uc(17) - LHestop  Stop the evolution when the helium luminosity exceeds this limit (if > 0)
real(double) :: uc(21)

! Optionally stop the code when the target mass is reached
logical :: stop_at_target_mass = .false.

end module stopping_conditions
