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
!!           Bit 4: (4) Indicates whether nucleosynthesis data is present
!!
!! ------------------------------------------------------------------------------
!<

subroutine remesh ( kh2, jch, bms, tm, p1, ecc, oa, jstar, jf )
   use real_kind
   use mesh
   use mesh_enc
   use extra_elements
   use constants
   use settings
   use atomic_data
   use test_variables
   use eostate_types
   use funcs1_interface
   use structure_variables
   use accretion_abundances
   use nucleosynthesis
   implicit none
   integer :: kh2, jch, jstar, jf
   real(double) :: bms, tm, p1, ecc, oa

   integer :: i,ik,ih
   real(double) :: ex_max,che,vma,qe,hpc,vd,q1,q2,vx,dk,qk,si
   real(double) :: new_h(nvar,NM), new_dh(nvar,NM)
   real(double) :: new_hnuc(50,NM)
   type(eostate) :: eos


   !Common blocks:
   real(double) :: var(nvar), dvar(nvar), fn1(nfunc)
   
   real(double) :: xh, xhe, xc, xn, xo, xne, xmg, xsi, xfe
   

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
   h(nMg24, 1:kh) = xmg
   if (use_mg24_eqn) then
      do ik=1, kh
         h(nmg24,ik) = max(0.0d0, 1.0 - (h(5,ik) + h(9, ik) + h(10, ik) + h(3, ik) + h(11, ik) + h(16, ik) + xfe + xsi))
      end do
   end if
   if ( jch >= 4 ) then
      h( 5,1:kh) = xh/vma
      h( 3,1:kh) = xo/vma
      h( 9,1:kh) = xhe/vma
      h(10,1:kh) = xc/vma
      h(11,1:kh) = xne/vma
      h(16,1:kh) = xn/vma
   end if
   h(nSi28, 1:kh) = xsi
   h(nFe56, 1:kh) = xfe

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

   qe = 0.0d0              ! Modified mesh spacing function
   mc(jstar) = tm          ! New mass after remesh
   var(:) = h(:, kh)
   dvar(:) = 0.0d0
   call funcs1 ( kh, -2, var(:), dvar(:), fn1(:), eos)
   hpc = dsqrt(eos%p/(cg*eos%rho*eos%rho))
   mc(jstar) = 3.5d-33*eos%rho*hpc**3
   vd = tm/h(4, 1)
   h(4, 1:kh) = vd*h(4, 1:kh)                ! Scale mass
   h(6, 1:kh) = 1.0d0                        ! Gradient of mesh spacing
   if ( iand(jf, 2) /= 2 ) h(12, 1:kh) = 1.0d0        ! Moment of Inertia (later)
   !H(13, 1:KH) = P1                          ! New rotation period
   h(17, 1:kh) = oa                          ! New orbital angular momentum
   h(18, 1:kh) = ecc                         ! New eccentricity
   h(20, 1:kh) = bms                         ! New total binary mass

   !     Rotation rate/rotational frequency
   h(13, 1:kh) = 2.0*cpi/(p1 * csday)
   dh(13, 1:kh) = -2.0*cpi/(p1**2 * csday) * dh(13, 1:kh)

   !     AGB MSF terms
   ct(11) = 0.1d0*eos%p
   ct(13) = 3.0d19
   ct(15) = 0.3d19

   ex_max = eos%ex
   if (ex_max < 1.0d-3) enc_parachute = h(8, 1) / h(4, 1)
   do ik = 1, kh
      var(:) = h(:, ik)
      call funcs1 ( ik, -2, var(:), dvar(:), fn1(:), eos )                 ! Calculate stuff
      qa(ik) = qq                            ! Store mesh spacing function
      ! Preliminary integration for M.I. and potential
      if ( .not. ( ik == 1 .or. iand(jf, 2)==2 ) ) then
         h(12, ik) = h(12, ik - 1) + r*r*m3/(dabs(qm))
         h(14, ik) = h(14, ik - 1) + phim/dabs(qm)
      end if
      ex_max = max(eos%ex, 0.0d0)
      ! Integration for L-dependent mesh spacing.
      !qe = qe + ct(2)*eos%ex*ht(4, ik)*1.0d-8
      !qa(ik) = qa(ik) + qe                   ! update meshspacing
   end do
   ! Find values of mesh spacing function at the external points
   ! Needed to calculate the new mesh, where the meshspacing gradient is
   ! constant.
   q1 = qa(1)
   q2 = (kh2 - 1.0d0)/(qa(kh) - qa(1))
   if ( jch >= 3 ) then
      ! If required, redistribute the mesh-points, either by changing
      ! their number or by changing the mesh-distribution function
      do ik = 1, kh
         vx = h(5, ik) + 1.0d-10             ! Fudge hydrogen fraction
         h(5, ik) = dlog(vx)
         qa(ik) = (qa(ik) - q1)*q2 + 1.0d0   ! Adjust meshspacing
      end do
      ih = 1
      do ik = 1, kh2
         dk = 0.0d0
         if ( ik == kh2 ) ih = kh
         if ( ik /= 1 .and. ik /= kh2 ) then
            ! Find the proper value for the meshspacing function at
            ! this meshpoint
            do i = 1, 50
               ! Sanity check: abort if we're running out of the mesh
               ! boundary
               if ( ih+1 > kh) then
                  write (0, *) 'remesh running outside mesh boundary, aborting'
                  stop
               end if
               if ( ik >= qa(ih + 1) ) ih = ih + 1
               if ( ik < qa(ih + 1) ) exit   ! Break loop
            end do
            dk = (ik - qa(ih))/(qa(ih + 1) - qa(ih))
         end if
         ! Linear interpolation for new H and DH
         if (ih < kh2) then
            new_h(:, ik) = h(:, ih) + dk*(h(:, ih + 1) - h(:, ih))
            new_dh(:, ik) = dh(:, ih) + dk*(dh(:, ih + 1) - dh(:, ih))
            if (nucleosynthesis_enabled) then
               new_hnuc(:, ik) = hnuc(Jstar, :, ih) + dk*(hnuc(Jstar, :, ih + 1) - hnuc(Jstar, :, ih))
            end if
         else
            new_h(:, ik) = h(:, ih)
            new_dh(:, ik) = dh(:, ih)
            if (nucleosynthesis_enabled) then
               new_hnuc(:, ik) = hnuc(Jstar, :, ih)
            end if
         endif
      end do
      q2 = q2 * (kh - 1.0) / (kh2 - 1.0)
      kh = kh2
      do ik = 1, kh
         ! Un-fudge the hydrogen abundance
         new_h(5, ik) = dexp(new_h(5, ik)) - 1.0d-10
         if ( new_h(5, ik) < 1.0d-5 ) new_h(5, ik) = 0.0d0
         h(:, ik) = new_h(:, ik)
         dh(:, ik) = new_dh(:, ik)
      end do
      if (nucleosynthesis_enabled) then
         if (kh_nuc < kh) then
            call allocate_nucleosynthesis_data(kh)
         end if
         hnuc(Jstar, :, 1:kh) = new_hnuc(:, 1:kh)
      end if
   end if

   qk = 1.0d0/q2                             ! Gradient of mesh spacing
   si = h(12, kh)                            ! Potential at the surface
   ! Some new variables that may not have been read in from stored model
   do ik = 1, kh
      h(6, ik) = qk                          ! Gradient of mesh spacing
      if (iand(jf,2) /= 2) then
         h(12, ik) = (si - h(12, ik))*dabs(qk)  ! Moment of inertia of interior
         h(14, ik) = - gmr - h(14, ik)*dabs(qk) ! Gravitational potential
         h(15, ik) = - gmr                      ! Potential at the stellar surface
      end if
   end do
   if (ex_max < 1.0d-3) enc_parachute = 1.1*h(8, 1) / h(4, 1)
   h(ntam, 1:kh) = h(12, 1:kh)*h(13, 1)

end subroutine remesh


