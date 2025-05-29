!> \file stars_structure_equations.f90   Contains equns1()
#include "assert.h"

#ifdef DEBUG
#define pure
#endif

pure subroutine equns1 ( jk, kl, kq, fn2, equv )
   use real_kind
   use mesh
   use mesh_enc
   use constants
   use settings
   use control
   use semi_implicit_variables
   use current_model_properties
   use accretion_abundances
   use step_functions
   use indices

   implicit none
   integer, intent(in) :: jk, kl, kq
   real(double), intent(in) :: fn2(3, NFUNC)
   real(double), intent(out) :: equv(NEQ)

   integer :: Jstar,ii
   real(double) :: dmu12,dmu23,dw12,dw23
   real(double) :: mt_smooth
   real(double) :: s12,s23
   real(double) :: si12,si23
   real(double) :: wta

   ! Difference-equation terms:
   real(double) :: bcp(3), bct(3), vp(3), vpk(3), vr(3), vrk(3)
   real(double) :: vt(3), vtk(3), L(3), lk(3), lq(3), mt(3), vm(3), vmk(3), sg(3)
   real(double) :: x1(3), x1t(3), x16(3), x16t(3), x4(3), x4t(3), x12(3)
   real(double) :: x12t(3), x20(3), x20t(3), bcm(3), vi(3), vik(3), phi(3), phik(3)
   real(double) :: bcf(3), bcs(3), bcph(3), x14(3), x14t(3), avmu(3), sgth(3), omega(3)
   real(double) :: omegat(3), sgam(3), si(3), adam(3), sam(3), samt(3)
   real(double) :: bca(3), bce(3), xi(3), xik(3), dlrk(3), bcmb(3), x24(3), x24t(3)
   real(double) :: x28(3), x28t(3), x56(3), x56t(3)
   real(double) :: am(3), amk(3)
   real(double) :: Vmc2(3), FU2(3), FV(3)
   real(double) :: fh(3), fhe(3), fc(3), fn(3), fo(3), fne(3), fmg(3), fsi(3), ffe(3)
   real(double) :: bcmp(3)
   real(double) :: macc

   ! equations and extra equations
   real(double) :: equ(neq)

   ! VAR(3),(2),(1): values at current, previous and anteprevious meshpoints.
   ! (3) is nearest centre if KL=0, KQ=1; nearest surface if KL=1, KQ=-1

   do Jstar = 1, ktw
   ! Copy input variables, star 1 or 2
   bcp(:)   = fn2(:, fn_idx_for_star(FN_BCP, Jstar))
   bct(:)   = fn2(:, fn_idx_for_star(FN_BCT, Jstar))
   vp(:)    = fn2(:, fn_idx_for_star(FN_VP, Jstar))
   vpk(:)   = fn2(:, fn_idx_for_star(FN_VPK, Jstar))
   vr(:)    = fn2(:, fn_idx_for_star(FN_VR, Jstar))
   vrk(:)   = fn2(:, fn_idx_for_star(FN_VRK, Jstar))
   vt(:)    = fn2(:, fn_idx_for_star(FN_VT, Jstar))
   vtk(:)   = fn2(:, fn_idx_for_star(FN_VTK, Jstar))
   L(:)     = fn2(:, fn_idx_for_star(FN_VL, Jstar))
   lk(:)    = fn2(:, fn_idx_for_star(FN_LK, Jstar))
   lq(:)    = fn2(:, fn_idx_for_star(FN_LQ, Jstar))
   mt(:)    = fn2(:, fn_idx_for_star(FN_MT, Jstar))
   vm(:)    = fn2(:, fn_idx_for_star(FN_VM, Jstar))
   vmk(:)   = fn2(:, fn_idx_for_star(FN_VMK, Jstar))
   sg(:)    = fn2(:, fn_idx_for_star(FN_SG, Jstar))
   x1(:)    = fn2(:, fn_idx_for_star(FN_X1, Jstar))
   x1t(:)   = fn2(:, fn_idx_for_star(FN_X1T, Jstar))
   x16(:)   = fn2(:, fn_idx_for_star(FN_X16, Jstar))
   x16t(:)  = fn2(:, fn_idx_for_star(FN_X16T, Jstar))
   x4(:)    = fn2(:, fn_idx_for_star(FN_X4, Jstar))
   x4t(:)   = fn2(:, fn_idx_for_star(FN_X4T, Jstar))
   x12(:)   = fn2(:, fn_idx_for_star(FN_X12, Jstar))
   x12t(:)  = fn2(:, fn_idx_for_star(FN_X12T, Jstar))
   x20(:)   = fn2(:, fn_idx_for_star(FN_X20, Jstar))
   x20t(:)  = fn2(:, fn_idx_for_star(FN_X20T, Jstar))
   bcm(:)   = fn2(:, fn_idx_for_star(FN_BCM, Jstar))
   vi(:)    = fn2(:, fn_idx_for_star(FN_VI, Jstar))
   vik(:)   = fn2(:, fn_idx_for_star(FN_VIK, Jstar))
   phi(:)   = fn2(:, fn_idx_for_star(FN_PHI, Jstar))
   phik(:)  = fn2(:, fn_idx_for_star(FN_PHIK, Jstar))
   bcf(:)   = fn2(:, fn_idx_for_star(FN_BCF, Jstar))
   bcs(:)   = fn2(:, fn_idx_for_star(FN_BCS, Jstar))
   bcph(:)  = fn2(:, fn_idx_for_star(FN_BCPH, Jstar))
   x14(:)   = fn2(:, fn_idx_for_star(FN_X14, Jstar))
   x14t(:)  = fn2(:, fn_idx_for_star(FN_X14T, Jstar))
   avmu(:)  = fn2(:, fn_idx_for_star(FN_AVMU, Jstar))
   sgth(:)  = fn2(:, fn_idx_for_star(FN_SGTH, Jstar))
   omega(:) = fn2(:, fn_idx_for_star(FN_OMEGA, Jstar))
   omegat(:)= fn2(:, fn_idx_for_star(FN_OMEGAT, Jstar))
   sam(:)   = fn2(:, fn_idx_for_star(FN_SAM, Jstar))
   samt(:)  = fn2(:, fn_idx_for_star(FN_SAMT, Jstar))
   adam(:)  = fn2(:, fn_idx_for_star(FN_ADAM, Jstar))
   sgam(:)  = fn2(:, fn_idx_for_star(FN_SGAM, Jstar))
   si(:)    = fn2(:, fn_idx_for_star(FN_SI, Jstar))

   am(:)    = fn2(:, fn_idx_for_star(fn_am, Jstar))
   amk(:)   = fn2(:, fn_idx_for_star(fn_amk, Jstar))

   Vmc2(:)  = fn2(:, fn_idx_for_star(fn_Vmc2, Jstar))
   FU2(:)   = fn2(:, fn_idx_for_star(fn_FU2k, Jstar))
   FV(:)    = fn2(:, fn_idx_for_star(fn_FV, Jstar))

   x24(:)   = fn2(:, fn_idx_for_star(fn_x24, Jstar))
   x24t(:)  = fn2(:, fn_idx_for_star(fn_x24t, Jstar))
   x28(:)   = fn2(:, fn_idx_for_star(fn_x28, Jstar))
   x28t(:)  = fn2(:, fn_idx_for_star(fn_x28t, Jstar))
   x56(:)   = fn2(:, fn_idx_for_star(fn_x56, Jstar))
   x56t(:)  = fn2(:, fn_idx_for_star(fn_x56t, Jstar))

   macc     = fn2(3, fn_idx_for_star(fn_macc, Jstar))

   fh       = fn2(3, fn_idx_for_star(fn_fh, Jstar))
   fhe      = fn2(3, fn_idx_for_star(fn_fhe, Jstar))
   fc       = fn2(3, fn_idx_for_star(fn_fc, Jstar))
   fn       = fn2(3, fn_idx_for_star(fn_fn, Jstar))
   fo       = fn2(3, fn_idx_for_star(fn_fo, Jstar))
   fne      = fn2(3, fn_idx_for_star(fn_fne, Jstar))
   fmg      = fn2(3, fn_idx_for_star(fn_fmg, Jstar))
   fsi      = fn2(3, fn_idx_for_star(fn_fsi, Jstar))
   ffe      = fn2(3, fn_idx_for_star(fn_ffe, Jstar))

   ! Input variables, binary orbit
   bca(:)   = fn2(:, fn_idx_for_star(index_orbit_fn_start + 1, Jstar))
   bce(:)   = fn2(:, fn_idx_for_star(index_orbit_fn_start + 2, Jstar))
   xi(:)    = fn2(:, fn_idx_for_star(index_orbit_fn_start + 3, Jstar))
   xik(:)   = fn2(:, fn_idx_for_star(index_orbit_fn_start + 4, Jstar))
   dlrk(:)  = fn2(:, fn_idx_for_star(index_orbit_fn_start + 5, Jstar))
   bcmb(:)  = fn2(:, fn_idx_for_star(index_orbit_fn_start + 6, Jstar))
   bcmp(:)  = fn2(:, fn_idx_for_star(FN_PMASS, Jstar))

   ! second-order difference equations at interior points
   if ( 3 <= jk + kl .and. jk + kl <= kh ) then

     ! Molecular weight gradient, for thermohaline mixing and angular momentum transport
     dmu12 = kq * (avmu(1)-avmu(2))
     dmu23 = kq * (avmu(2)-avmu(3))

     ! Gradient in rotational velocity, for rotational shear
     dw12 = (omega(2) - omega(1))**2
     dw23 = (omega(3) - omega(2))**2

     ! Combined diffusion coefficients for chemical mixing
     ! Convection and thermohaline mixing
     mt_smooth = 0.0d0
     s12 = 0.5d0*(sg(1) + sg(2)) - mixing_fudge*pstv(kq*mt(2), mt_smooth) &
          + 0.5d0*(sgth(1)+sgth(2))*pstv(dmu12, 0.0d0)
     s23 = 0.5d0*(sg(2) + sg(3)) - mixing_fudge*pstv(-kq*mt(3), mt_smooth)&
          + 0.5d0*(sgth(2)+sgth(3))*pstv(dmu23, 0.0d0)

     equ(EQN_H1)   = s23*(x1(3)  - x1(2))  - s12*(x1(2) - x1(1))  - x1t(2)
     equ(EQN_HE4)  = s23*(x4(3)  - x4(2))  - s12*(x4(2) - x4(1))  - x4t(2)
     equ(EQN_C12)  = s23*(x12(3) - x12(2)) - s12*(x12(2) -x12(1)) - x12t(2)
     equ(EQN_N14)  = s23*(x14(3) - x14(2)) - s12*(x14(2) -x14(1)) - x14t(2)
     equ(EQN_O16)  = s23*(x16(3) - x16(2)) - s12*(x16(2) -x16(1)) - x16t(2)
     equ(EQN_NE20) = s23*(x20(3) - x20(2)) - s12*(x20(2) -x20(1)) - x20t(2)
     equ(EQN_MG24) = s23*(x24(3) - x24(2)) - s12*(x24(2) -x24(1)) - x24t(2)
     equ(EQN_SI28) = s23*(x28(3) - x28(2)) - s12*(x28(2) -x28(1)) - x28t(2)
     equ(EQN_FE56) = s23*(x56(3) - x56(2)) - s12*(x56(2) -x56(1)) - x56t(2)

     ! Add advection terms (contributions from gravitational settling)
     if (cgrs > 0.0d0) then
        equ(EQN_H1)   = equ(EQN_H1)   - kq*(fh(3)  - fh(2))
        equ(EQN_HE4)  = equ(EQN_HE4)  - kq*(fhe(3) - fhe(2))
        equ(EQN_C12)  = equ(EQN_C12)  - kq*(fc(3)  - fc(2))
        equ(EQN_N14)  = equ(EQN_N14)  - kq*(fn(3)  - fn(2))
        equ(EQN_O16)  = equ(EQN_O16)  - kq*(fo(3)  - fo(2))
        equ(EQN_NE20) = equ(EQN_NE20) - kq*(fne(3) - fne(2))
        equ(EQN_MG24) = equ(EQN_MG24) - kq*(fmg(3) - fmg(2))
        equ(EQN_SI28) = equ(EQN_SI28) - kq*(fsi(3) - fsi(2))
        equ(EQN_FE56) = equ(EQN_FE56) - kq*(ffe(3) - ffe(2))
     end if

     ! Advection-diffusion equation for angular momentum
     !equ(EQN_OMEGA) = s23*(omega(3) - omega(2)) - s12*(omega(2) - omega(1))
     s12  = 0.5d0*(sgam(1) + sgam(2)) + 0.5d0*(sgth(1)+sgth(2))*pstv(dmu12, 0.0d0)
     s23  = 0.5d0*(sgam(2) + sgam(3)) + 0.5d0*(sgth(2)+sgth(3))*pstv(dmu23, 0.0d0)
     si12 = 0.5d0*(sgam(1)*si(1) + sgam(2)*si(2)) + 0.5d0*(sgth(1)*si(1)+sgth(2)*si(2))*pstv(dmu12, 0.0d0)
     si23 = 0.5d0*(sgam(2)*si(2) + sgam(3)*si(3)) + 0.5d0*(sgth(2)*si(2)+sgth(3)*si(3))*pstv(dmu23, 0.0d0)
     s12  = 0.0d0
     s23  = 0.0d0
     equ(EQN_OMEGA) = -max(kq*adam(3), 0.0d0) + max(-kq*adam(1), 0.0d0) - adam(2) & ! Advection
                      + s23*(sam(3) - sam(2)) - s12*(sam(2) - sam(1))&              ! Diffusion
                      + si23*(omega(3) - omega(2)) - si12*(omega(2) - omega(1))&    ! Shear
                      + max(kq*mt(2), 0.0d0)*(sam(2)-sam(1)) - max(-kq*mt(3), 0.0d0)*(sam(3)-sam(2))& ! Advection due moving mesh
                      - samt(2)
   end if

   ! Next-to-surface boundary conditions for second-order equations
   if ( jk + kl == 2 ) then
     ! Modified equations for accreting of matter with different composition

     ! Thermohaline mixing of surface layer
     if (kl == 0) then    ! (3) is nearest to centre, so we have the grid points | 3 | 2 | AC |
        s12 = 0.5d0*(sg(1) + sg(2)) - mixing_fudge*pstv(kq*mt(2), 0.0d0)
        s23 = 0.5d0*(sg(2) + sg(3)) - mixing_fudge*pstv(-kq*mt(3), 0.0d0)

        equ(EQN_H1)   = s23*(x1(3)  - x1(2))  - s12*(x1(2)  - xac(1,Jstar)) - x1t(2)
        equ(EQN_HE4)  = s23*(x4(3)  - x4(2))  - s12*(x4(2)  - xac(2,Jstar)) - x4t(2)
        equ(EQN_C12)  = s23*(x12(3) - x12(2)) - s12*(x12(2) - xac(3,Jstar)) - x12t(2)
        equ(EQN_N14)  = s23*(x14(3) - x14(2)) - s12*(x14(2) - xac(4,Jstar)) - x14t(2)
        equ(EQN_O16)  = s23*(x16(3) - x16(2)) - s12*(x16(2) - xac(5,Jstar)) - x16t(2)
        equ(EQN_NE20) = s23*(x20(3) - x20(2)) - s12*(x20(2) - xac(6,Jstar)) - x20t(2)
        equ(EQN_MG24) = s23*(x24(3) - x24(2)) - s12*(x24(2) - xac(7,Jstar)) - x24t(2)
        equ(EQN_SI28) = s23*(x28(3) - x28(2)) - s12*(x28(2) - xac(8,Jstar)) - x28t(2)
        equ(EQN_FE56) = s23*(x56(3) - x56(2)) - s12*(x56(2) - xac(9,Jstar)) - x56t(2)
     else                 ! (3) is nearest to surface, we have the grid points | 1 | 2 | 3 | AC |
        s12 = 0.5d0*(sg(2) + sg(3)) - pstv(kq*mt(3), 0.0d0)&
             + 0.5d0*(sgth(2)+sgth(3))*pstv(kq * (avmu(2)-avmu(3)), 0.0d0)
        s23 = -macc

        equ(EQN_H1)   = s23*(xac(1,Jstar) - x1(3))  - s12*(x1(3)  - x1(2))  - x1t(3)
        equ(EQN_HE4)  = s23*(xac(2,Jstar) - x4(3))  - s12*(x4(3)  - x4(2))  - x4t(3)
        equ(EQN_C12)  = s23*(xac(3,Jstar) - x12(3)) - s12*(x12(3) - x12(2)) - x12t(3)
        equ(EQN_N14)  = s23*(xac(4,Jstar) - x14(3)) - s12*(x14(3) - x14(2)) - x14t(3)
        equ(EQN_O16)  = s23*(xac(5,Jstar) - x16(3)) - s12*(x16(3) - x16(2)) - x16t(3)
        equ(EQN_NE20) = s23*(xac(6,Jstar) - x20(3)) - s12*(x20(3) - x20(2)) - x20t(3)
        equ(EQN_MG24) = s23*(xac(7,Jstar) - x24(3)) - s12*(x24(3) - x24(2)) - x24t(3)
        equ(EQN_SI28) = s23*(xac(8,Jstar) - x28(3)) - s12*(x28(3) - x28(2)) - x28t(3)
        equ(EQN_FE56) = s23*(xac(9,Jstar) - x56(3)) - s12*(x56(3) - x56(2)) - x56t(3)
     end if

     ! Advection terms (from gravitational settling)
     if (cgrs > 0.0d0) then
        equ(EQN_H1)   = equ(EQN_H1)   - kq*(fh(3)  - fh(2))
        equ(EQN_HE4)  = equ(EQN_HE4)  - kq*(fhe(3) - fhe(2))
        equ(EQN_C12)  = equ(EQN_C12)  - kq*(fc(3)  - fc(2))
        equ(EQN_N14)  = equ(EQN_N14)  - kq*(fn(3)  - fn(2))
        equ(EQN_O16)  = equ(EQN_O16)  - kq*(fo(3)  - fo(2))
        equ(EQN_NE20) = equ(EQN_NE20) - kq*(fne(3) - fne(2))
        equ(EQN_MG24) = equ(EQN_MG24) - kq*(fmg(3) - fmg(2))
        equ(EQN_SI28) = equ(EQN_SI28) - kq*(fsi(3) - fsi(2))
        equ(EQN_FE56) = equ(EQN_FE56) - kq*(ffe(3) - ffe(2))
     end if

     ! Angular momentum transport boundary condition - rigid rotation
     s23  = sgam(2)       + sgth(2)*pstv(dmu23, 0.0d0)
     si23 = sgam(2)*si(2) + sgth(2)*si(2)*pstv(dmu23, 0.0d0)
     s23  = 0.0d0
     equ(EQN_OMEGA) = -max(kq*adam(3), 0.0d0) + max(-kq*adam(1), 0.0d0) - adam(2) & ! Advection
                      + s23*(sam(3) - sam(2))            &                          ! Diffusion
                      + si23*(omega(3) - omega(2))       &                          ! Shear
                      - max(-mt(3), 0.0d0)*sam(3)        &                          ! Mass loss
                      - samt(2)
   end if

   ! Next-to-central boundary conditions for second-order equations
   if ( jk + kl == kh + 1 ) then
     s23 = kq*0.5d0*(sg(2)+sg(3)) + kq*0.5d0*(sgth(2)+sgth(3))*pstv(kq * (avmu(2)-avmu(3)), 0.0d0)
     equ(EQN_H1)   = s23*(x1(3)  - x1(2))  + x1t(3 - kl)
     equ(EQN_HE4)  = s23*(x4(3)  - x4(2))  + x4t(3 - kl)
     equ(EQN_C12)  = s23*(x12(3) - x12(2)) + x12t(3 - kl)
     equ(EQN_N14)  = s23*(x14(3) - x14(2)) + x14t(3 - kl)
     equ(EQN_O16)  = s23*(x16(3) - x16(2)) + x16t(3 - kl)
     equ(EQN_NE20) = s23*(x20(3) - x20(2)) + x20t(3 - kl)
     equ(EQN_MG24) = s23*(x24(3) - x24(2)) + x24t(3 - kl)
     equ(EQN_SI28) = s23*(x28(3) - x28(2)) + x28t(3 - kl)
     equ(EQN_FE56) = s23*(x56(3) - x56(2)) + x56t(3 - kl)

     ! Angular momentum transport
     !s23 = kq*0.5d0*(sgam(2)*si(2) + sgam(3)*si(3))
     !equ(EQN_OMEGA) = s23*(omega(3) - omega(2)) + omegat(3 - kl)
     !si23 = 0.5d0*(sgam(2)*si(2) + sgam(3)*si(3))
     !s23 = 0.0d0
     !equ(EQN_OMEGA) = - max(-kq*mt(3), 0.0d0)*(sam(3)-sam(2)) - samt(3 - KL)
     !equ(EQN_OMEGA) = - samt(3-KL)
     !equ(EQN_OMEGA) = omega(3) - omega(2)
     !equ(EQN_OMEGA) = si23*(omega(3) - omega(2))

     !equ(EQN_OMEGA) = adam(3) - adam(2) &                                           ! Advection
     !                 + s23*(sam(3) - sam(2)) &                                     ! Diffusion
     !                 + si23*(omega(3) - omega(2)) &                                ! Shear
     !                 + samt(3-kl)

     ! Use solid body rotation; since the angular momentum of the inner gridpoint is 0 we need to be
     ! careful with it: it could well act as an infinite sink for angular momentum
     equ(EQN_OMEGA) = omega(3) - omega(2)
   end if

   ! The sum of all abundances should be conserved (1, in fact)
   equ(EQN_SUMX) = equ(EQN_H1)  + equ(EQN_HE4)  + equ(EQN_C12)  + equ(EQN_N14) +  &
                   equ(EQN_O16) + equ(EQN_NE20) + equ(EQN_MG24) + equ(EQN_SI28) + equ(EQN_FE56)

   ! First-order difference equations at interior points.
   ! KL=0, KQ= 1 means that XX(3) is deeper in star than XX(2);
   ! KL=1, KQ=-1 means that XX(3) nearer to surface than XX(2).
   ! XXK is sort-of dXX/dK, except that its sign must be flipped if KL=1,
   ! so a factor of KQ is put in.
   if ( 2.le.jk .and. jk.le.kh ) then
     ii = 2*Jstar - 3
     wta = 0.5d0*kq
     equ(EQN_PRES)   = vp(3)  - vp(2)  - wta*(vpk(3)  + vpk(2))   ! Pressure
     equ(EQN_RADIUS) = vr(3)  - vr(2)  - wta*(vrk(3)  + vrk(2))   ! Radius
     equ(EQN_TEMP)   = vt(3)  - vt(2)  - wta*(vtk(3)  + vtk(2))   ! Temperature
     equ(EQN_MASS)   = vm(3)  - vm(2)  - wta*(vmk(3)  + vmk(2))   ! Mass
     equ(EQN_INERT)  = vi(3)  - vi(2)  - wta*(vik(3)  + vik(2))   ! Moment of inertia
     equ(EQN_PHI)    = phi(3) - phi(2) - wta*(phik(3) + phik(2))  ! Gravitational potential
     equ(EQN_XI)     = xi(3)  - xi(2)  - wta*(xik(3)  + xik(2))   ! Mass-transfer rate
     equ(EQN_TAM)    = am(3)  - am(2)  - wta*(amk(3)  + amk(2))   ! Total angular momentum
     equ(EQN_LUM)    = L(3)   - L(2)   - wta*(lk(3) + lk(2)) &    ! Luminosity
                     - lq(2)*pstv(kq*mt(2),0.0d0) + lq(3)*pstv(-kq*mt(3),0.0d0) &     ! Advection term
                     + ii*kq*dlrk(2 + kl)                                             ! Heat trabsfer in contact
   end if

   ! Surface boundary conditions for first-order equations and `eigenvalues'
   ! For the numbering of equ(), see the section 'The boundary conditions' in the manual
   if ( jk == 1 ) then
     equ(SBC_MDOT)  = bcm(3)   ! Mass loss
     equ(SBC_PRES)  = bcp(3)   ! Pressure
     equ(SBC_TEMP)  = bct(3)   ! Temperature
     equ(SBC_PHI)   = bcf(3)   ! Potential
     equ(SBC_SSPIN) = bcs(3)   ! Spin, solid body
     equ(SBC_PHIS)  = bcph(3)  ! Surface potential (eigenvalue)
     equ(SBC_OAM)   = bca(3)   ! Orbital angular momentum
     equ(SBC_ECC)   = bce(3)   ! Orbital eccentricity
     equ(SBC_BMASS) = bcmb(3)  ! Binary mass
     equ(SBC_PMASS) = bcmp(3)  ! Primary mass
   end if

   ! Central boundary conditions for first-order equations
   ! For the numbering of equ(), see the section 'The boundary conditions' in the manual
   if ( jk == kh + 1 ) then
     equ(CBC_MASS)   = vm(3)   ! Mass
     equ(CBC_LUM)    = L(3)    ! Luminosity
     equ(CBC_RADIUS) = vr(3)   ! Radius
     equ(CBC_INERT)  = vi(3)   ! Moment of inertia
     equ(CBC_XI)     = xi(3)   ! Mass-transfer rate
     equ(CBC_TAM)    = am(3)   ! Total angular momentum
   end if

   ! Copy results to output

   ! Equations, star 1 or 2
   forall (ii = 1:nestar) equv(idx_for_star(ii, Jstar)) = equ(ii)

   ! Equations, binary orbit
   forall (ii = INDEX_ORBIT_EQN_START:INDEX_ORBIT_EQN_START + NVBIN) equv(ii) = equ(ii)

   end do  ! do Jstar = 1, ktw
end subroutine equns1

