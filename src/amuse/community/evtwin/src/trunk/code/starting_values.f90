
!> module starting_values0
!!
!! Contains variables from the former common block t0
!< 

module starting_values0
   use real_kind
   implicit none
   
   real(double) :: t0sm     ! Primary mass
   real(double) :: t0dty    ! Time step
   real(double) :: t0age    ! Age
   real(double) :: t0per    ! Orbital period
   real(double) :: t0bms    ! Binary mass
   real(double) :: t0ecc    ! Orbital eccentricity
   real(double) :: t0p      ! Rotational period
   real(double) :: t0enc    ! Energy-generation term
   
end module starting_values0

