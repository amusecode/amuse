module atomic_data
   use real_kind
   
   implicit none
   ! Constants in partition function for H2 from Vardya (1960), Webbink (1975)
   real(double) :: ch2(4)
   real(double) :: chi(26,9)     ! Ionisation potentials [eV]
   real(double) :: com(27)       ! Statistical weights
   real(double) :: can(9)        ! Atomic mass number (in units of AMU)
   real(double) :: cbn(9)        ! Baryon number (protons+neutrons, != can)
   integer :: kzn(9)       ! Number of protons/electrons

   ! These are for computational convenience
   real(double) :: dkzn(9)       ! Number of protons/electrons
   real(double) :: lcan(9)       ! log(can**1.5), used a few times in the EoS
end module

