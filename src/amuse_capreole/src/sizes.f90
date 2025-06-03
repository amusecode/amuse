!>
!! \brief This module contains basic size parameter definitions
!!
!! Module for Capreole
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2010-02-04
!!
!! This module is also accepted by the F compiler (Dec 9, 2003)

module sizes

  private

  integer,public,parameter :: nrOfDim=3        !< number of dimensions
  integer,public,parameter :: mbc=2            !< number of ghost cells
  integer,public,parameter :: neuler=2+NrOfDim !< number of Euler equations
  integer,public,parameter :: neq=neuler+2     !< number of equations (Euler + advected quantities)
  ! Indices of state array
  integer,public,parameter :: RHO=1 !< Indices of state array: density
  integer,public,parameter :: RHVX=2 !< Indices of state array: x momentum density
  integer,public,parameter :: RHVY=3 !< Indices of state array: y momentum density
  integer,public,parameter :: RHVZ=4 !< Indices of state array: z momentum density
  integer,public,parameter :: EN=5 !< Indices of state array: energy density
  !integer,public,parameter :: TRACER1=6
  integer,public,parameter :: XHI=6 !< Indices of state array: neutral hydrogen fraction
  integer,public,parameter :: XHII=7 !< Indices of state array: ionized hydrogen fraction

  ! The following define constants for identifying the coordinate system
  integer,public,parameter :: CART =1 !< constant for defining cartesian coordinate system
  integer,public,parameter :: CYL  =2 !< constant for defining cylindrical coordinate system
  integer,public,parameter :: SPHER=3 !< constant for defining spherical coordinate system

end module sizes
