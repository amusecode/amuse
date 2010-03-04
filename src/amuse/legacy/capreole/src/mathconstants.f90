!>
!! \brief This module contains mathematical constants
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! Author: Garrelt Mellema
!!
!! Date: 2003-12-09
!!
!! This module is also accepted by the F compiler (Dec 9, 2003)
!<
module mathconstants

  use precision, only: dp

  implicit none

  private

  !> the number pi
  real(kind=dp),public,parameter :: pi=3.141592654

end module mathconstants
