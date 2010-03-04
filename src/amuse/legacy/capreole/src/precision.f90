!>
!!    \brief This module contains the floating point precision definitions 
!!    in a machine
!!    independent way
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! Author: Garrelt Mellema
!!
!! Date: 2003-12-09
!!
!! This module is also accepted by the F compiler (Dec 9, 2003)
!<
module precision

  ! Module for Capreole / C2-Ray (f90)
  ! Author: Garrelt Mellema
  ! Date: 2003-012-09
  ! This module is also accepted by the F compiler (Dec 9, 2003)

  ! This module contains the precision definitions in a machine
  ! independent way

  implicit none

  private

  !> This is single precision (4B)
  integer,public,parameter :: si=selected_real_kind(6,37) 
  !> This is double precision (8B)
  integer,public,parameter :: dp=selected_real_kind(15,307) 
  !> This is integer (4B)
  integer,public,parameter :: li=selected_int_kind(9)
  real(kind=dp),public,parameter :: tiny_dp=tiny(1.0_dp) !< smallest dp

end module precision
