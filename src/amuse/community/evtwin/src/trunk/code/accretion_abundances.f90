!> module accretion_abundances:
!! 
!! Abundances of accreted material for star 1 and star 2 in a binary.
!! Also used during post-He flash construction to preserve the composition profile.
!< 

module accretion_abundances
   use real_kind
   implicit none
   
   ! Default accretion abundances (from init.dat)
   real(double) :: x1ac      ! Hydrogen abundance
   real(double) :: x4ac      ! Helium abundance
   real(double) :: x12ac     ! Carbon abundance
   real(double) :: x14ac     ! Neon abundance
   real(double) :: x16ac     ! Oxygen abundance
   real(double) :: x20ac     ! Neon abundance
   real(double) :: x24ac     ! Magnesium abundance
   real(double) :: x28ac     ! Silicon abundance
   real(double) :: x56ac     ! Iron abundance
   
   ! Actual abundances used for accretion onto star 1 or star 2
   ! May be different from the above if the code is keeping track of the true abundances during accretion.
   real(double) :: xac(9,2)  ! Stores the above elements, for star 1, 2
   
end module accretion_abundances

