!>
!! \brief This module contains atomic constants and parameter definitions
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2010-02-04 (2003-06-01)
!!
!! \b Version: cgs
!!
!! This module is also accepted by the F compiler (Dec 9, 2003)
!<

module atomic

  use precision, only: dp
  use scaling, only: scvelo
  use cgsconstants, only: kb, m_p

  private

  !> Boltzmann constant over proton mass (gas constant)
  real(kind=dp),public,parameter :: boltzm = kb/(m_p*scvelo*scvelo)

  !> adiabatic index
  real(kind=dp),public,parameter :: gamma = 5.0_dp/3.0_dp
  real(kind=dp),public,parameter :: gamma1 = gamma - 1.0_dp

end module atomic
