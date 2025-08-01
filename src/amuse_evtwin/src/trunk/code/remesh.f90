!> ------------------------------------------------------------------------------
!! REMESH
!!  Adjust the mesh spacing of a model.
!! ------------------------------------------------------------------------------
!!  Input:
!!     KH2 - new number of meshpoints
!!           KH in the common block holds the current number of meshpoints
!!     JCH - How to remesh:
!!           1, 2: Only initialise new variables (the options are equivalent)
!!           3:    Calculate new mesh spacing and redistribute meshpoints
!!           4:    As (3), but also homogenise composition profile
!!     BMS - New total binary mass            [Msun]
!!      TM - New primary mass                 [Msun]
!!      P1 - New rotational period for star   [days]
!!     ECC - New orbital eccentricity
!!      OA - New orbital angular momentum     [g cm^2/s]
!!   JSTAR - Which star in the binary to remesh (1 or 2)
!!      JF - Flag indicating which variables were present in the input file
!!           and which need to be recalculated (bitfield):
!!           Bit 1: (1) Unused
!!           Bit 2: (2) if set, take potential and moment of inertia from model
!!           Bit 3: (4) Indicates whether DH(:,:) is present in the input file
!!           Bit 4: (8) Indicates whether nucleosynthesis data is present
!!
!! ------------------------------------------------------------------------------
#include "assert.h"

subroutine remesh ( kh2, jch, bms, tm, p1, ecc, oa, jstar, jf )
   use real_kind
   use mesh
   use mesh_enc
   use constants
   use settings
   use atomic_data
   use test_variables
   use eostate_types
   use structure_functions
   use structure_variables
   use accretion_abundances
   use nucleosynthesis
   use interpolate
   use indices
   implicit none
   integer, intent(in) :: kh2, jch, jstar, jf
   real(double), intent(in) :: bms, tm, p1, ecc, oa
   real(double) :: r

   integer :: i,ik,ikk
   real(double) :: ex_max,che,vma,qe,hpc,vd,q1,q2,qk
   type(eostate) :: eos


   !Common blocks:
   real(double) :: var(nvar), dvar(nvar), fn1(nfunc)
   real(double) :: vars(nvar), varc(nvar)
   real(double) :: qa(NM)
   type(interpolate_t) :: interpol_h(NVSTAR + NVBIN)
   type(interpolate_t) :: interpol_dh(NVSTAR + NVBIN)
   type(interpolate_t) :: interpol_hnuc(50)

   real(double) :: xh, xhe, xc, xn, xo, xne, xmg, xsi, xfe
   real(double) :: qq            ! Determines mesh-point metric: mesh-point interval
   real(double) :: qm            ! Determines mesh-point metric: derivative of mesh spacing function wrt mass
   real(double) :: phim          ! Derivative of gravitational potential with respect to m**(2/3)
   real(double) :: gmr           ! Effective gravity at the surface(?)
   real(double) :: m3            ! m^(1/3), m is the mass coordinate
   real(double) :: mc1           ! Core mass for star 1, backup
   real(double) :: nmf

   ! Sanity checks
   assert(kh2 <= max_nm)
   assert(kh <= max_nm)
   assert(allocated(h))
   assert(allocated(sx))
   assert(oa >= 0.0d0)
   assert(ecc >= 0.0d0)
   assert(ecc <= 1.0d0)
   assert(p1 > 0.0d0)
   assert(bms > tm)

   ! Record initial age
   age0 = age

   ! Set initial composition.
   ! The composition variables are NOT the actual mass fractions if we
   ! use non-integer atomic masses, so we have to compute what they are
   ! The composition variables used in the code are baryon number fractions
   ! We calculate this even if we don't want to do a ZAMS run because we need to
   ! know the baryon number densities of Fe, Si and Mg.
   che = 1.0d0 - ch - czs
   cn = 1.0d0 - cc - co - cne - cmg - csi - cfe
   xh = ch*cbn(1)/can(1)
   xhe = che*cbn(2)/can(2)
   xc = cc*czs*cbn(3)/can(3)
   xn = cn*czs*cbn(4)/can(4)
   xo = co*czs*cbn(5)/can(5)
   xne = cne*czs*cbn(6)/can(6)
   xmg = cmg*czs
   xsi = csi*max(czs, 1.0d-4)
   xfe = cfe*max(czs, 1.0d-4)
   vma = xh + xhe + xc + xn + xo + xne + xmg + xsi + xfe
   xfe = xfe / vma
   xsi = xsi / vma
   xmg = xmg / vma

   ! Mg24 is never stored in the input model!
   xmg = max(1.0 - (xh+xhe+xc+xn+xo+xne)/vma-(xsi+xfe), 0.0d0)
   h(VAR_MG24, 1:kh) = xmg
   if (use_mg24_eqn) then
      do ik=1, kh
         xmg = h(VAR_H1,ik) + h(VAR_HE4, ik) + h(VAR_C12, ik) + h(VAR_O16, ik) + h(VAR_NE20, ik) + h(VAR_N14, ik) + xfe + xsi
         h(VAR_MG24,ik) = max(0.0d0, 1.0 - xmg)
      end do
   end if
   if ( jch >= 4 ) then
      h(VAR_H1,  1:kh) = xh/vma
      h(VAR_O16, 1:kh) = xo/vma
      h(VAR_HE4, 1:kh) = xhe/vma
      h(VAR_C12, 1:kh) = xc/vma
      h(VAR_NE20,1:kh) = xne/vma
      h(VAR_N14, 1:kh) = xn/vma
   end if
   h(VAR_SI28, 1:kh) = xsi
   h(VAR_FE56, 1:kh) = xfe

   ! Now we must also convert the abundances of the accreted material. If not
   ! set from init.dat, set from initial abundances.
   if (x1ac < 0.0) x1ac = ch*cbn(1)/can(1)
   if (x4ac < 0.0) x4ac = che*cbn(2)/can(2)
   if (x12ac < 0.0) x12ac = cc*czs*cbn(3)/can(3)
   if (x14ac < 0.0) x14ac = cn*czs*cbn(4)/can(4)
   if (x16ac < 0.0) x16ac = co*czs*cbn(5)/can(5)
   if (x20ac < 0.0) x20ac = cne*czs*cbn(6)/can(6)
   if (x24ac < 0.0) x24ac = cmg*czs*cbn(7)/can(7)
   if (x28ac < 0.0) x28ac = csi*czs*cbn(8)/can(8)
   if (x56ac < 0.0) x56ac = cfe*czs*cbn(9)/can(9)
   ! make sure XH is 1-everything else and abundancies sum to 1
   x1ac = max(0.0d0, 1.0d0 -(x4ac+x12ac+x14ac+x16ac+x20ac+x24ac+csi*czs+cfe*czs))

   xh = x1ac*cbn(1)/can(1)
   xhe = x4ac*cbn(2)/can(2)
   xc = x12ac*cbn(3)/can(3)
   xn = x14ac*cbn(4)/can(4)
   xo = x16ac*cbn(5)/can(5)
   xne = x20ac*cbn(6)/can(6)
   xmg = x24ac*cbn(7)/can(7)
   xsi = x28ac*cbn(8)/can(8)
   xfe = x56ac*cbn(9)/can(9)

   vma = xh + xhe + xc + xn + xo + xne + xmg + xsi + xfe

   x1ac  = xh / vma
   x4ac  = xhe / vma
   x12ac = xc / vma
   x14ac = xn / vma
   x16ac = xo / vma
   x20ac = xne / vma
   x24ac = xmg / vma
   x28ac = xsi / vma
   x56ac = xfe / vma

   ! Initialise accretion abundances for both stars
   xac(1, 1:2) = x1ac
   xac(2, 1:2) = x4ac
   xac(3, 1:2) = x12ac
   xac(4, 1:2) = x14ac
   xac(5, 1:2) = x16ac
   xac(6, 1:2) = x20ac
   xac(7, 1:2) = x24ac
   xac(8, 1:2) = x28ac
   xac(9, 1:2) = x56ac

   ! Backup core mass for star 1
   mc1 = mc(1)

   ! Get central properties of the star
   qe = 0.0d0              ! Modified mesh spacing function
   mc(jstar) = tm          ! Initialise core mass (it needs to have a non-zero value)
   var(:) = h(:, kh)
   dvar(:) = 0.0d0
   call funcs1 ( kh, -2, var(:), dvar(:), fn1(:), eos, px=sx(:, 2))
   hpc = sqrt(eos%p/(cg*eos%rho*eos%rho))
   mc(jstar) = 3.5d-33*eos%rho*hpc**3

   ! Scale the mass of the star
   vd = tm/h(VAR_MASS, 1)
   h(VAR_MASS, 1:kh) = vd*h(VAR_MASS, 1:kh)

   h(VAR_QK, 1:kh) = 1.0d0                            ! Gradient of mesh spacing
   if ( iand(jf, 2) /= 2 ) h(VAR_INERT, 1:kh) = 1.0d0 ! Moment of Inertia (later)
   if ( iand(jf, 2) /= 2 ) h(VAR_PHI, 1:kh)   = 1.0d0 ! Potential
   !H(13, 1:KH) = P1                                  ! New rotation period
   h(VAR_HORB, 1:kh) = oa                             ! New orbital angular momentum
   h(VAR_ECC, 1:kh) = ecc                             ! New eccentricity
   h(VAR_BMASS, 1:kh) = bms                           ! New total binary mass

   !     Rotation rate/rotational frequency
   h(VAR_OMEGA, 1:kh) = 2.0*cpi/(p1 * csday)
   dh(VAR_OMEGA, 1:kh) = -2.0*cpi/(p1**2 * csday) * dh(VAR_OMEGA, 1:kh)

   !     AGB MSF terms
   if (ct(11) == 0.0d0) ct(11) = 0.1d0*eos%p
   if (ct(13) == 0.0d0) ct(13) = 3.0d19
   if (ct(15) == 0.0d0) ct(15) = 0.3d19

   ! Set core mass as core mass for star 1.
   ! This is a bit of a hack that is needed because funcs1 will always think that the current star is star 1
   mc(1) = mc(jstar)

   ex_max = eos%ex
   do ik = 1, kh
      ikk = kh + 2 - ik
      var(:) = h(:, ik)
      call funcs1 ( ik, -2, var(:), dvar(:), fn1(:), eos, px=sx(:,ikk)) ! Calculate stuff
      qq = sx(85, ikk)
      qm = sx(86, ikk)
      phim = sx(87, ikk)
      m3 = sx(89, ikk)
      qa(ik) = qq                            ! Store mesh spacing function
      r   = sqrt(abs(exp(2.0d0*var(7)) - ct(8)))
      ! Preliminary integration for M.I. and potential
      if ( .not. ( ik == 1 .or. iand(jf, 2)==2 ) ) then
         h(VAR_PHI,   ik) = h(VAR_PHI,   ik - 1) + phim/abs(qm)
      end if
      ex_max = max(eos%ex, 0.0d0)
      ! Integration for L-dependent mesh spacing.
      !qe = qe + ct(2)*eos%ex*ht(4, ik)*1.0d-8
      !qa(ik) = qa(ik) + qe                   ! update meshspacing
   end do

   ! Store gravitational potential
   gmr = sx(88, kh+1)

   ! Restore core mass of star 1 if this was really star 2 */
   if (jstar == 2) mc(1) = mc1

   ! Scale factor in number of meshpoints
   nmf = real(kh2 / kh)

   ! Find values of mesh spacing function at the external points
   ! Needed to calculate the new mesh, where the meshspacing gradient is
   ! constant.
   qk = (qa(kh) - qa(1)) / (kh - 1.0d0)
   if ( jch >= 3 ) then
      ! Normalise mesh spacing function
      q1 = qa(1)
      qa(1) = 1.0d0
      do ik = 2, kh
         qa(ik) = (qa(ik) - q1) + 1.0d0
      end do
      q2 = qa(kh)

      ! Create interpolation splines for all variables
      do i = 1, NVSTAR+NVBIN
         call make_interpolation_table(kh, qa(1:kh), h(i, 1:kh), interpol_h(i))
         call make_interpolation_table(kh, qa(1:kh), dh(i, 1:kh), interpol_dh(i))
      end do
      if (nucleosynthesis_enabled) then
         do i = 1, 50
            call make_interpolation_table(kh, qa(1:kh), hnuc(Jstar, i, 1:kh), interpol_hnuc(i))
         end do
      end if
      ! Safe the centre and surface properties (these should not be affected by remeshing)
      vars(:) = h(:, 1)
      varc(:) = h(:, kh)

      ! Calculate new mesh spacing function
      qa(1) = 1.0d0
      qa(kh2) = q2
      qk = qk * (kh - 1.0d0) / (kh2 - 1.0d0)
      kh = kh2
      do ik = 2, kh-1
         qa(ik) = qa(ik-1) + qk
      end do

      ! Calculate stellar properties for the new mesh
      do ik = 1, kh
         do i = 1, NVSTAR+NVBIN
            h(i, ik) = evaluate_interpolation_table(qa(ik), interpol_h(i))
            dh(i, ik) = evaluate_interpolation_table(qa(ik), interpol_dh(i))
         end do
         if (nucleosynthesis_enabled) then
            do i = 1, 50
               hnuc(Jstar, i, ik) = evaluate_interpolation_table(qa(ik), interpol_hnuc(i))
            end do
         end if
      end do
      ! Restore centre and surface properties
      h(:, 1) = vars(:)
      h(:, kh) = varc(:)
   end if

   ! Some new variables that may not have been read in from stored model
   h(VAR_QK, 1:kh) = qk                                                    ! Gradient of mesh spacing
   if ( iand(jf, 2)/=2 ) then
      h(VAR_PHIS,  1:kh) = -gmr                                            ! Potential at the stellar surface
      forall (ik = 1:kh) h(VAR_PHI,   ik) = -gmr - h(VAR_PHI, ik)*abs(qk)*nmf  ! Gravitational potential

      ! Integrate moment of inertia, the straightforward way
      h(VAR_INERT, kh) = 0.0d0
      do ik = kh-1, 1, -1
         var(:) = h(:, ik)
         call funcs1 ( ik, -2, var(:), dvar(:), fn1(:), eos, px=sx(:,ikk))
         h(VAR_INERT, ik) = h(VAR_INERT, ik + 1) - fn1(fn_vik)
      end do
   end if

   h(VAR_TAM, 1:kh) = h(VAR_INERT, 1:kh)*h(VAR_OMEGA, 1)
   h(VAR_PMASS, 1:kh) = h(VAR_MASS, 1)

end subroutine remesh


