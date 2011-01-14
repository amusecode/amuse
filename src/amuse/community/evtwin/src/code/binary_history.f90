!> binary_history:
!! 
!! This module contains the history of the primary and the orbit when doing binary evolution
!! in non-TWIN mode
!! 
!! \todo Give variables somewhat longer names for easier recognition and to avoid conflicts (?)
!< 

module binary_history
   use real_kind
   use mesh
   
   implicit none
   real(double) :: hpr(NVAR,NM)  ! H() of primary ???
   real(double) :: ms(9999)      ! Mass of ... ???
   
   real(double) :: sect(9999)    ! Time ???
   real(double) :: sdt(9999)     ! dt ???
   real(double) :: scm(9999)     ! Companion mass (of which star?) ???
   real(double) :: sang(9999)    ! Orbital AM ???
   real(double) :: se(9999)      ! Orbital excentricity ???
   
   real(double) :: ww(16)        ! Dummy ???
   real(double) :: wx(15)        ! Dummy ???
   
end module binary_history

