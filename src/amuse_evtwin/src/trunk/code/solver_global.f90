! -*- Mode: Fortran90; tab-width: 3 -*-
module solver_global_variables
   use real_kind
   use nucleosynthesis, only: nvar_nuc
   use mesh
   use indices, only: neq
   implicit none

   integer :: block_vars = 0
   integer :: block_nmesh = 0

   real(double), allocatable :: S(:, :)
   real(double), allocatable :: C(:,:,:)
   real(double), allocatable :: xd(:)
   real(double), allocatable :: diffx(:,:)
   real(double), allocatable :: er(:)          ! Typical value of variables
   real(double)              :: es(NEQ) = 0.0d0! Typical values, previous iter
   real(double), allocatable :: eqn_scale(:)   ! Typical value of eqns
   real(double), allocatable :: sbc_scale(:)   ! Typical value of sbcs
   real(double), allocatable :: ddh(:,:)
   !  real(double), allocatable :: GRADF(:)

   ! Default scalings for variables. Normally 1, but these can be
   ! overwritten by other parts of the code (eg. printb)
   real(double), allocatable :: default_er(:)

   ! Quantities needed for "linesearch" algorithm, see equations (9.7.7)
   ! through (9.7.9) in Numerical Recipes (2nd edition)
   real(double) :: residue         ! called g in NR
   real(double) :: prev_residue    ! Value of RESIDUE on previous iter
   real(double) :: residue_deriv   ! called g' in NR

   ! How many equations can be solved for; init.dat only allows input of
   ! up to 40, but this can be extended in principle. For
   ! nucleosynthesis, we need at least 50 (45 actually).
   integer, parameter :: nkd = 50
   integer :: keq, kvb, kvc
   integer :: kq, kee
   integer :: kj2, kj5, kj6, kj10, kj12, ki4
   integer :: ke1, ke2, ke3, ke4, kbc, kev, kl, jh1, jh2, jh3
   integer :: ki1, ki2, ki3, kj1, kj11, kj3, kj4, kj7, kj8, kj9
   integer :: ik_first, ik_last
   integer :: kd(3*nkd)    ! 3 sets of NKD integers, for NKD variables

   real(double), allocatable :: func(:,:)!(NFUNC, NM)
   real(double), allocatable :: dfunc(:,:,:)!(NFUNC, NVAR, NM)
   
   integer :: solver_output_unit(3) = (/ 1, 35, 36 /)

end module solver_global_variables


