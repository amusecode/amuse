module model_initialiser

contains

   ! remesher:
   ! Function based on the origianl REMESH function
   ! Input (in COMMON blocks):
   !  H(,) - array of independent variables (in unnamed COMMON)
   !  DH(,) - list of changes in independent variables (in unnamed COMMON)
   !  KH - current number of meshpoints in H(,)
   ! Input options:
   !  KH2 - new number of meshpoints for this model
   !  JCH - switch to determine whether to construct new mesh spacing function
   !        and initialise composition or not.
   !  BMS - Total mass in the binary (EV)
   !  TM  - New mass of the star
   !  P1  - New rotational period of the star (solid body)
   !  ECC - New eccentricity of the binary
   !  OA  - New orbital angular momentum
   !  JSTAR - Labels which if the stars to initislise variables for
   !  JF - Switch to decide which variables to recompute
   ! TODO: this could be split up into different subroutines for each of the
   ! individual tasks.
   subroutine remesher( kh2, jch, bms, tm, p1, ecc, oa, jstar, jf )
      use real_kind
      use constants
      use mesh
      use init_dat
      use settings
      use control
      use atomic_data
      use test_variables
      use eostate_types
      use structure_functions
      use current_model_properties
      use structure_variables
      use accretion_abundances
      use indices
      
      implicit none
      integer, intent(in) :: kh2, jch, jstar, jf
      real(double), intent(in) :: bms, tm, p1, ecc, oa
      integer :: nm_current, nm_next, nm_target
      integer :: ik, ikk, ih, i
      integer :: jo
      integer :: ksv, kt5
      type(init_dat_settings) :: initdat
      real(double) :: nh(nvar,nm), ndh(nvar,nm), nndh(nvar,nm)
      real(double) :: q1, q2, dk, dty, vd, dtb, ageb, pcrit, wcrit, w1
      real(double) :: si, vma, hpc
      logical :: equilibrium
      logical :: newmesh
      real(double) :: var(nvar), dvar(nvar), fn1(nfunc)
      real(double) :: qa(NM)
      type(eostate) :: eos

      real(double) :: xh0, xhe0, xc0, xn0, xo0, xne0, xmg0, xsi0, xfe0
      real(double) :: che
      real(double) :: xh, xhe, xc, xn, xo, xne, xmg, xsi, xfe
      real(double) :: r

      real(double) :: qq            ! Determines mesh-point metric: mesh-point interval
      real(double) :: qm            ! Determines mesh-point metric: derivative of mesh spacing function wrt mass
      real(double) :: phim          ! Derivative of gravitational potential with respect to m**(2/3)
      real(double) :: gmr           ! Effective gravity at the surface(?)
      real(double) :: m3            ! m^(1/3), m is the mass coordinate


      ! Backup current settings so we can restore them when we're done
      dtb = dt
      ageb = age
      call push_init_dat(initdat, kh2, ksv, kt5, jch)
      ! Change settings to a reasonable set of defaults
      call load_basic_init_dat(ik, ksv, kt5, ik)
      kop = initdat%kop
      kx = 0; ky = 0; kz = 0; kth = 0
      cmi = 0.0
      crd = 0.0d0
      kt1 = 100
      kt2 = 0
      kt3 = 0
      kt4 = 100
      kt5 = 100
      ksv = 100
      joc = 1
      jter = 0

      ! Set initial composition.
      ! The composition variables are NOT the actual mass fractions if we
      ! use non-integer atomic masses, so we have to compute what they are
      ! The composition variables used in the code are baryon number fractions
      ! We calculate this even if we don't want to do a ZAMS run because we need to
      ! know the baryon number densities of Fe, Si and Mg.
      ! Note that the actual abundances of Si and Fe are forced to by non-zero.
      ! This is because the EoS becomes poorly defined if there are not enough free
      ! electrons. It has no impact on the opacity and only a small impact on the
      ! mean molecular weight and the mass loss rate.
      che = 1.0d0 - ch - czs
      cn = 1.0d0 - cc - co - cne - cmg - csi - cfe

      xh0 = ch*cbn(1)/can(1)
      xhe0 = che*cbn(2)/can(2)
      xc0 = cc*czs*cbn(3)/can(3)
      xn0 = cn*czs*cbn(4)/can(4)
      xo0 = co*czs*cbn(5)/can(5)
      xne0 = cne*czs*cbn(6)/can(6)
      xmg0 = cmg*czs*cbn(7)/can(7)
      xsi0 = csi*czs*cbn(8)/can(8)
      xfe0 = cfe*czs*cbn(9)/can(9)

      xh = xh0
      xhe = xhe0
      xc = xc0
      xn = xn0
      xo = xo0
      xne = xne0
      xmg = xmg0
      xsi = max(xsi0, csi*1.0d-4)
      xfe = max(xfe0, cfe*1.0d-4)

      vma = xh + xhe + xc + xn + xo + xne + xmg + xsi + xfe
      xh = xh / vma
      xhe = xhe / vma
      xc = xc / vma
      xn = xn / vma
      xo = xo / vma
      xne = xne / vma
      xfe = xfe / vma
      xsi = xsi / vma
      xmg = xmg / vma

      ! Initialise composition variables
      h(VAR_MG24, 1:kh) = xmg
      h(VAR_SI28, 1:kh) = xsi
      h(VAR_FE56, 1:kh) = xfe
      if ( jch >= 4 ) then
         h(VAR_H1,1:kh) = xh
         h(VAR_O16,1:kh) = xo
         h(VAR_HE4,1:kh) = xhe
         h(VAR_C12,1:kh) = xc
         h(VAR_NE20,1:kh) = xne
         h(VAR_N14,1:kh) = xn
      end if
      ! We should always do this for Mg24, since that's never stored
      if (use_mg24_eqn) then
         do ik=1, kh
            xmg = h(VAR_H1,ik) + h(VAR_HE4, ik) + h(VAR_C12, ik) + h(VAR_O16, ik) + h(VAR_NE20, ik) + h(VAR_N14, ik) + xfe + xsi
            h(VAR_MG24,ik) = max(0.0d0, 1.0 - xmg)
         end do
      end if

      ! Now we must also convert the abundances of the accreted material. If not
      ! set from init.dat, set from initial abundances.
      !> \todo FIXME: this will cause problems if we ever need to call REMESH twice in
      !! the same run
      !<
      x1ac = x1ac*cbn(1)/can(1)
      x4ac = x4ac*cbn(2)/can(2)
      x12ac = x12ac*cbn(3)/can(3)
      x14ac = x14ac*cbn(4)/can(4)
      x16ac = x16ac*cbn(5)/can(5)
      x20ac = x20ac*cbn(6)/can(6)
      x24ac = x24ac*cbn(7)/can(7)
      if (x1ac < 0.0) x1ac = xh0
      if (x4ac < 0.0) x4ac = xhe0
      if (x12ac < 0.0) x12ac = xc0
      if (x14ac < 0.0) x14ac = xn0
      if (x16ac < 0.0) x16ac = xo0
      if (x20ac < 0.0) x20ac = xne0
      if (x24ac < 0.0) x24ac = xmg0
      vma = x1ac + x4ac + x12ac + x14ac + x16ac + x20ac + x24ac + xfe + xsi

      x1ac  = x1ac  / vma
      x4ac  = x4ac  / vma
      x12ac = x12ac / vma
      x14ac = x14ac / vma
      x16ac = x16ac / vma
      x20ac = x20ac / vma
      x24ac = x24ac / vma
      ! make sure XH is 1-everything else and abundancies sum to 1
      x1ac = max(0.0d0, 1.0d0-(x4ac+x12ac+x14ac+x16ac+x20ac+x24ac+xfe0+xsi0))

      ! Initialise accretion abundances for both stars
      xac(1, 1:2) = x1ac
      xac(2, 1:2) = x4ac
      xac(3, 1:2) = x12ac
      xac(4, 1:2) = x14ac
      xac(5, 1:2) = x16ac
      xac(6, 1:2) = x20ac
      xac(7, 1:2) = x24ac

      ! Set initial values of some other variables
      ! Typical mass-scale for the interior (needed for mesh spacing function)
      mc(jstar) = tm          ! New mass after remesh
      var(:) = h(:, kh)
      dvar(:) = 0.0d0
      call funcs1 ( kh, -2, var(:), dvar(:), fn1(:), eos, px=sx(:,2))
      hpc = sqrt(eos%p/(cg * eos%rho * eos%rho))
      mc(jstar) = 3.5d-33 * eos%rho * hpc**3

      ! Initialise binary (orbital) parameters
      h(VAR_HORB, 1:kh) = oa                             ! New orbital angular momentum
      h(VAR_ECC, 1:kh) = ecc                             ! New eccentricity
      h(VAR_BMASS, 1:kh) = bms                           ! New total binary mass

      ! First: change the number of meshpoints or the mesh spacing function
      ! Determine if we need to calculate a new mesh (independent of the number
      ! of meshpoints)
      newmesh = .false.
      if (jch>3) newmesh = .true.
      ! Set new number of meshpoints
      nm_current = kh
      nm_target = kh2
      nm_next = nm_target
      print *, 'nremesh from', nm_current, 'to', nm_target
      do while(newmesh .or. nm_current /= nm_next)
         newmesh = .false.
         print *, 'trying ', nm_next
         ! Store old model, so we can go back if needed
         nh(:,1:nm_current) = h(:,1:nm_current)
         ndh(:,1:nm_current) = dh(:,1:nm_current)
         ! Find values of mesh spacing function
         do ik=1, nm_current
            var(:) = h(:, ik)
            call funcs1 ( ik, -2, var(:), dvar(:), fn1(:))                 ! Calculate stuff
            qa(ik) = qq                            ! Store mesh spacing function
         end do
         ! Interpolate model onto new mesh
         ! Find values of mesh spacing function at the external points
         ! Needed to calculate the new mesh, where the meshspacing gradient is
         ! constant.
         q1 = qa(1)
         q2 = (nm_next - 1.0d0)/(qa(nm_current) - qa(1))
         do ik = 1, nm_current
            qa(ik) = (qa(ik) - q1)*q2 + 1.0d0   ! Adjust meshspacing
         end do
         ih = 1
         do ik = 1, nm_next
            dk = 0.0d0
            if ( ik == nm_next ) ih = nm_current
            if ( ik /= 1 .and. ik /= nm_next ) then
               ! Find the proper value for the meshspacing function at
               ! this meshpoint
               do i = 1, 50
                  ! Sanity check: abort if we're running out of the mesh
                  ! boundary
                  if ( ih+1 > nm_current) then
                     write (0, *) &
                          'remesh running outside mesh boundary, aborting'
                     stop
                  end if
                  if ( ik >= qa(ih + 1) ) ih = ih + 1
                  if ( ik < qa(ih + 1) ) exit   ! Break loop
               end do
               dk = (ik - qa(ih))/(qa(ih + 1) - qa(ih))
            end if
            ! Linear interpolation for new H and DH
            h(:, ik) = nh(:, ih) + dk*(nh(:, ih + 1) - nh(:, ih))
            nndh(:, ik) = ndh(:, ih) + dk*(ndh(:, ih + 1) - ndh(:, ih))
         end do
         !H(6, 1:KH) = 1.0/Q2              ! Gradient of mesh spacing

         ! Now see if the model will converge properly if we let the code iterate
         dty = dt/csy
         jo = 0
         jnn = 0
         kh = nm_next
         call printb ( jo, 1, 22 )
         age = age - dty
         call nextdt ( dty, jo, 22 )
         jnn = 1
         dh(:, 1:nm_next) = 0.0d0
         call solver(20, id, kt5, jo)
         if (jo == 0) then
            print *, 'converged ok'
            ! If yes, pick next number of meshpoints
            nm_current = nm_next
            nm_next = nm_target
            h(:,1:nm_current) = h(:,1:nm_current) + dh(:,1:nm_current)
         else
            ! If no, pick a smaller number of meshpoints in between the current value
            ! and the target and try again.
            nm_next = (nm_current+nm_next)/2
            print *, 'cannot converge, reduce to ', nm_next
            ! Restore backup copies of H and DH
            h(:,1:nm_current) = nh(:,1:nm_current)
            dh(:,1:nm_current) = ndh(:,1:nm_current)
         end if
      end do
      print *, 'nremesh finished with ', nm_current, '(wanted ', nm_target,')'
      if (nm_current < nm_target) then
         print *, '*** nremesh failed ***'
         stop
      end if
      dh(:, 1:nm_current) = nndh(:,1:nm_current)
      dh(:, 1:nm_current) = 0.0
      kh = nm_current

      ! Second: scale the mass
      vd = tm/h(VAR_MASS, 1)
      vd = 1.0
      print *, 'scaling mass by factor', vd
      do while (abs(1.0d0 - vd) > 0.1)
         vd = max(0.9d0,min(vd, 1.1d0))
         h(VAR_MASS, 1:kh) = vd*h(VAR_MASS, 1:kh)       ! Scale mass

         kth = 1
         jhold = 4
         equilibrium = equilibrate_model(kt5)
         if (.not. equilibrium) then !H(:,1:nm_current) = NH(:,1:nm_current)
            print *, '*** failed ***'
            stop
         end if

         vd = tm/h(VAR_MASS, 1)
      end do
      h(VAR_MASS, 1:kh) = vd*h(VAR_MASS, 1:kh)          ! Scale mass


      ! Third: scale the surface rotation rate
      ! Make the whole star rotate with the surface rate, if desired
      ! Convert rotational period -> rotation rate
      forall (ik=1:kh) h(VAR_OMEGA, ik) = 2.0*cpi/(h(VAR_OMEGA, ik) * csday)
      if (start_with_rigid_rotation) h(VAR_OMEGA,2:kh) = h(VAR_OMEGA,1)
      ! Now scale the surface rate, similar to the way the mass is scaled
      ! If we forced the star to rigid rotation before then this will set the
      ! rotation profile throughout the entire star
      wcrit = sqrt(cg*tm/exp(3*h(VAR_LNR, 1)))
      pcrit = 2.0*cpi/wcrit/csday
      w1 = 2.0*cpi/p1/csday
      print *, 'scaling rotation rate by factor', w1/h(13, 1)
      ! For rotation rates within 1/4 of critical be a bit more careful. Here we
      ! need to approach the desired rate smoothly, adjusting the structure of
      ! the star at each step.
      !> \todo FIXME: This should be more akin to the bit that updates the code on the
      !! new mesh, ie, not necessarily bring the star into equilibrium.
      !<
      if (wcrit/w1 < 4.0) then
         call set_solid_rotation
         print *, 'Period', p1, 'close to critical rate', pcrit
         vd = 0.25*wcrit/h(VAR_OMEGA,1)
         h(VAR_OMEGA, 1:kh) = vd*h(VAR_OMEGA, 1:kh)     ! Scale rotation rate
         kth = 1
         equilibrium = converge_model(kt5)
         do while (equilibrium .and. abs(w1 - h(VAR_OMEGA,1)) > 1.0d-6)
            vd = max(0.9d0,min(vd, 1.1d0))
            h(VAR_OMEGA, 1:kh) = vd*h(VAR_OMEGA, 1:kh)
            jhold = 4
            equilibrium = converge_model(kt5)
            vd = w1/h(VAR_OMEGA, 1)
            if (.not. equilibrium) then
               print *, '*** failed ***'
               stop
            end if
         end do
      end if
      vd = w1/h(VAR_OMEGA, 1)
      h(VAR_OMEGA, 1:kh) = vd*h(VAR_OMEGA, 1:kh)        ! Scale rotational period

      ! Compute moment of inertia and surface potential
      if (jf /= 2) then
         q2 = h(VAR_QK,1)
         h(VAR_QK, 1:kh) = 1.0d0
         h(VAR_INERT, 1:kh) = 1.0d0              ! Moment of Inertia
         h(VAR_PHI, 1:kh) = 0.0d0              ! Gravitational potential
         do ik = 2, kh
            ikk = kh + 2 - ik
            var(:) = h(:, ik)
            call funcs1 ( ik, -2, var(:), dvar(:), fn1(:), px=sx(:,ikk))
            qq = sx(85, ikk)
            qm = sx(86, ikk)
            phim = sx(87, ikk)
            m3 = sx(89, ikk)

            r = sqrt(abs(exp(2.0d0*var(7)) - ct(8)))
            h(VAR_INERT, ik) = h(VAR_INERT, ik - 1) + r*r*m3/(abs(qm))
            h(VAR_PHI, ik) = h(VAR_PHI, ik - 1) + phim/abs(qm)
         end do
         gmr = sx(88, kh+1)
         h(VAR_QK, 1:kh) = q2
         h(VAR_PHIS, 1:kh) = - gmr              ! Potential at the stellar surface
         si = h(VAR_INERT, kh)                   ! Total moment of inertia

         do ik = 1, kh
            h(VAR_INERT, ik) = (si - h(VAR_INERT, ik))*abs(q2)  ! Moment of inertia of interior
            h(VAR_PHI, ik) = - gmr - h(VAR_PHI, ik)*abs(q2) ! Gravitational potential
         end do
      end if
      ! Total angular momentum integration
      h(VAR_TAM, 1:kh) = h(VAR_INERT, 1:kh)*h(VAR_OMEGA, 1)

      if (relax_loaded_model) then
         kth = 1
         jhold = 4
         print *, 'equilibrating...'
         equilibrium = equilibrate_model(kt5)
         if (.not. equilibrium) then
            print *, '*** failed ***'
            stop
         end if
         print *, 'done'
         dh(:, 1:kh) = 0.0
      end if

      ! Restore old init.dat
      call pop_init_dat(initdat, ik, ksv, kt5, ik)
      dt = dtb
      age = ageb

      print *, 'Remesh done'
   end subroutine remesher



   function converge_model(kt5)
      use real_kind
      use mesh
      use constants
      use test_variables
      
      implicit none
      logical :: converge_model
      integer, intent(in) :: kt5
      real(double) :: dty
      integer :: jo

      jo = 0
      dty = dt/csy
      call solver(20, id, kt5, jo)
      if (jo == 0) then
         age = 0.0d0
         call printb ( jo, 1, 22 )
         h(:,1:kh) = h(:,1:kh) + dh(:,1:kh)
         call nextdt ( dty, jo, 22 )
      end if
      converge_model = (jo == 0)

   end function converge_model



   function equilibrate_model(kt5)
      use real_kind
      use mesh
      use constants
      use test_variables
      
      implicit none
      logical :: equilibrate_model
      integer, intent(in) :: kt5
      real(double) :: dty
      integer :: i, jo

      jo = 0
      do i=1, 40
         dty = dt/csy
         call solver(20, id, kt5, jo)
         if (jo /= 0) exit
         age = 0.0d0
         call printb ( jo, 1, 22 )
         h(:,1:kh) = h(:,1:kh) + dh(:,1:kh)
         call nextdt ( dty, jo, 22 )
         if (abs(lth) < 1.0d-8 .or. lth < 0.0d0) exit
      end do
      equilibrate_model = .false.
      if (lth < 1.0d-6 .and. jo == 0) equilibrate_model = .true.

   end function equilibrate_model

end module model_initialiser


