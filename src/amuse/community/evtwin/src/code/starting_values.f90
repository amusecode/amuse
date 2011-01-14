
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




!> module starting_values1
!!
!! Contains variables from the former common block tn1
!!
!! \todo This module is only used in ms2bss() and read_init_run().
!! In the former routine, only t1p is used!
!< 

module starting_values1
   use real_kind
   implicit none
   
   real(double) :: t1sm     ! Primary mass
   real(double) :: t1dty    ! Time step
   real(double) :: t1age    ! Age
   real(double) :: t1per    ! Orbital period
   real(double) :: t1bms    ! Binary mass
   real(double) :: t1ecc    ! Orbital eccentricity
   real(double) :: t1p      ! Rotational period
   real(double) :: t1enc    ! Energy-generation term
   
end module starting_values1

