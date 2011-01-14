module semi_implicit_variables
   use real_kind
   use mesh
   implicit none
   
   ! Second order corrections for composition equations
   real(double) :: ddX1(NM)     ! Hydrogen
   real(double) :: ddX4(NM)     ! Helium
   real(double) :: ddX12(NM)    ! Carbon
   real(double) :: ddX14(NM)    ! Nitrogen
   real(double) :: ddX16(NM)    ! Oxygen
   real(double) :: ddX20(NM)    ! Neon
   real(double) :: ddX24(NM)    ! Magnesium
   
   ! Second order corrections for luminosity equation
   real(double) :: ddL(NM)
   
   ! Interpolation for Lagrangian entropy terms (duplicated from mesh_enc)
   real(double), save :: TM(NM)
   real(double), save :: TSa(NM)
   real(double), save :: TSb(NM)
   real(double), save :: TSc(NM)
   real(double), save :: TSd(NM)
   
   
contains
   
   
   !> ------------------------------------------------------------------------------
   !!  GET_SM
   !!   Get the interppolated value of the entropy S for the given mass
   !!   coordinate.
   !!
   !!   \todo FIXME: make this work for star 1 and for star 2.
   !!
   !! ------------------------------------------------------------------------------
   !!  Input:
   !!     M  - Mass coordinate.
   !!  Return value:
   !!     The value of S at this mass.
   !! ------------------------------------------------------------------------------
   !<
   function get_sm(m)
      use real_kind
      use mesh
      use interpolate
      
      implicit none
      real(double), intent(in) :: m
      real(double) :: get_sm
      
      get_sm = iptable_eval(kh, m, TM(:), TSa(:), TSb(:), TSc(:), TSd(:))
   end function get_sm
   
end module semi_implicit_variables



! ------------------------------------------------------------------------------
!  COMP_SEMI_IMPLICIT_QUANTITIES
!   computes the values of quantities that are updated in between
!   iterations of the solver (hence semi-implicit).
! ------------------------------------------------------------------------------
!  Input:
!     //          - H and DH
!     structure_variables     - SX
!  Output:
!     ddXn        - Second order corrections to diffusion equation of species n
! ------------------------------------------------------------------------------
subroutine comp_semi_implicit_quantities
   use real_kind
   use mesh
   use semi_implicit_variables
   use control
   use interpolate
   use structure_variables
   
   implicit none
   ! Functions:
   real(double), external :: fracstep
   
   ! Local variables:
   real(double) :: X1(KH), X4(KH), X12(KH), X14(KH), X16(KH), X20(KH), X24(KH)
   integer :: k, kk
   
   
   if (.not. apply_second_order_corrections) return
   
   ! Get current value for variables, H+DH
   X1(1:KH)  = H( 5, 1:KH) + DH( 5, 1:KH)
   X4(1:KH)  = H( 9, 1:KH) + DH( 9, 1:KH)
   X12(1:KH) = H(10, 1:KH) + DH(10, 1:KH)
   X14(1:KH) = H(16, 1:KH) + DH(16, 1:KH)
   X16(1:KH) = H( 3, 1:KH) + DH( 3, 1:KH)
   X20(1:KH) = H(11, 1:KH) + DH(11, 1:KH)
   X24(1:KH) = H(41, 1:KH) + DH(41, 1:KH)
   
   ! Clear corrections
   ddX1(1:KH+1)  = 0.0d0    ! Hydrogen
   ddX4(1:KH+1)  = 0.0d0    ! Helium
   ddX12(1:KH+1) = 0.0d0    ! Carbon
   ddX14(1:KH+1) = 0.0d0    ! Nitrogen
   ddX16(1:KH+1) = 0.0d0    ! Oxygen
   ddX20(1:KH+1) = 0.0d0    ! Neon
   ddX24(1:KH+1) = 0.0d0    ! Magnesium
   
   ! Calculate second order corrections at each meshpoint
   do k=2, kh-1
      ddX1(k)  = fracstep( ( X1(k+1)- X1(k))*( X1(k)- X1(k-1)),  X1(k+1) -X1(k-1) )
      ddX4(k)  = fracstep( ( X4(k+1)- X4(k))*( X4(k)- X4(k-1)),  X4(k+1)- X4(k-1) )
      ddX12(k) = fracstep( (X12(k+1)-X12(k))*(X12(k)-X12(k-1)), X12(k+1)-X12(k-1) )
      ddX14(k) = fracstep( (X14(k+1)-X14(k))*(X14(k)-X14(k-1)), X14(k+1)-X14(k-1) )
      ddX16(k) = fracstep( (X16(k+1)-X16(k))*(X16(k)-X16(k-1)), X16(k+1)-X16(k-1) )
      ddX20(k) = fracstep( (X20(k+1)-X20(k))*(X20(k)-X20(k-1)), X20(k+1)-X20(k-1) )
      ddX24(k) = fracstep( (X24(k+1)-X24(k))*(X24(k)-X24(k-1)), X24(k+1)-X24(k-1) )
   end do
   
   ! Advection term in luminosity equation
   ddL(1:KH+1) = 0.0d0
   do k=2, kh-1
      kk = kh + 2 - k
      ddL(k)  = fracstep( ( SX(58, kk+1)- SX(58, kk))*( SX(58, kk)- SX(58, kk-1)),  SX(58, kk+1) -SX(58, kk-1) )
   end do
   
   ! Interpolation functions for entropy profile
   !> \todo FIXME: does SX actually contain the values for the previous
   !! timestep? It may contain "trial values" for the current timestep
   !! (in fact it should for the ddL() terms above to make sense).
   !< 
   do k=1, kh
      kk = kh + 2 - k
      TM(k) = H(4, k)
      TSa(k) = SX(28, kk)
   end do
   
   call iptable_init (kh, TM(:), TSa(:), TSb(:), TSc(:), TSd(:))
   
end subroutine comp_semi_implicit_quantities



! ------------------------------------------------------------------------------
! FRACSTEP
! ------------------------------------------------------------------------------
! Input:
!  F - Numerator of fraction to be caculated
!  G - Denominator of fraction to be calculated
! Returns:
!  F/G if F>0, 0 otherwise
! ------------------------------------------------------------------------------
function fracstep (f, g)
   use real_kind
   
   implicit none
   real(double), intent(in) :: f, g
   real(double) :: fracstep
   
   fracstep = 0.0d0
   if (f > 0.0d0) fracstep = f / g
end function fracstep

