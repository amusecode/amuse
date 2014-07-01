module polytrope
   use real_kind

   integer, parameter :: NMP = 499!1000

   real(double), private :: avg_rho, rho_c, Z
   real(double), private :: a, k, n, PBC, MBC
   real(double), private :: rho(NMP), P(NMP), T(NMP), lnf(NMP)
   real(double), private :: kappa(NMP), GRADAp(NMP)
   real(double), private :: r(0:NMP), m(0:NMP), L(0:NMP), Inert(0:NMP)
   real(double), private :: Q(NMP), QQ(NMP), QC
   real(double), private :: Lnuc, TKH

contains
!  Construct a polytropic stellar model (n=3/2) with a given mass and radius.
!  Because the surface of the star needs to satisfy the photospheric boundary
!  condition the surface of the star does not coincide with the surface of the
!  polytrope. We adjust the radius of the polytrope until we match the
!  photospheric boudnary condition at the desired radius. We then adjust the
!  mass of the polytrope until the mass within the photosphere matches the
!  desired stelar mass.
!
!  dm and dr are intial guesses for the required corrections (dm may be 0
!  initially but dr should not be).
   subroutine get_polytrope(mass, radius, xa, dm, dr)
      use constants
      use settings
      use poly32
      implicit none
      character (512) :: str
      real(double), intent(in) :: mass, radius, xa(9)
      real(double), intent(inout) :: dm, dr
      real(double) :: ddm, ddr
      real(double) :: BCM(3), BCP(3), MM(2,2), invM(2,2), detM
      real(double) :: GRADR
      integer :: i, iter, IJ
      logical :: verbose = .false.

      ddr = 0.0d0
      ddm = 0.0d0

! Now adjust the mass and radius of the polytrope until we match the
! photospheric boundary when we truncate the polytrope at the desired
! stellar radius. The second condition is that we match the total stellar
! mass within the photoshpere.
! The entire polytrope will therefore be larger and more massive than the
! star itself.
      if (radius+dr < 0.0d0) dr = 1.0d-3*radius
      if (radius+dr < 0.0d0) stop
      do iter = 1, 100
         call compute_polytrope_constants(mass+dm, radius+dr)
         call setup_mesh_spacing(radius, radius+dr, xa)
         call integrate_polytrope(mass, xa)

         if (verbose) write(6, "(1X, 1P, I2, ' Mass and radius:', 2E12.5)") &
                 iter, m(NMP)/CMSN, r(NMP)/CRSN

!     Compute numerical derivatives
         ddm = 1.0d-6
         ddr = 1.0d-6
         BCM(1) = MBC
         BCP(1) = PBC

         call compute_polytrope_constants(mass+dm+ddm, radius+dr)
         call integrate_polytrope(mass, xa)
         BCM(2) = MBC
         BCP(2) = PBC

         call compute_polytrope_constants(mass+dm, radius+dr+ddr)
         call integrate_polytrope(mass, xa)
         BCM(3) = MBC
         BCP(3) = PBC

!        Newton-Raphson step: find corrections by multiplying the residuals
!        with the inverse of the matrix of derivatives
         if(verbose) write(6, "(1X, 1P, '   Deviations      ', 2E12.5)")  BCP(1), BCM(1)
         MM(1,1) = (BCP(2) - BCP(1)) / ddm; MM(1,2) = (BCP(3) - BCP(1)) / ddr
         MM(2,1) = (BCM(2) - BCM(1)) / ddm; MM(2,2) = (BCM(3) - BCM(1)) / ddr
         detM = MM(2,2)*MM(1,1) - MM(1,2)*MM(2,1)
         invM(1,1) =  MM(2,2) / detM
         invM(1,2) = -MM(1,2) / detM
         invM(2,2) =  MM(1,1) / detM
         invM(2,1) = -MM(2,1) / detM

         ddm = invM(1,1) * BCP(1) + invM(1,2) * BCM(1)
         ddr = invM(2,1) * BCP(1) + invM(2,2) * BCM(1)
         if(verbose) write(6, "(1X, 1P, '   Corrections     ', 2E12.5)")  ddm, ddr

         dm = dm - ddm
         dr = dr - ddr 

!        Break out if we've converged
         if ( (ddm**2+ddr**2) < 1.0d-14) exit
      end do
      call compute_polytrope_constants(mass+dm, radius+dr)
      call setup_mesh_spacing(radius, radius+dr, xa)
      call integrate_polytrope(mass, xa)

      if (verbose) then
         print *, 'Polytrope mass and radius:', (mass+dm)/CMSN, (radius+dr)/CRSN
         print *, 'Stellar mass and radius:  ', m(NMP)/CMSN, r(NMP)/CRSN
         print *, 'Luminosity (total, nuc.): ', L(NMP)/CLSN, Lnuc/CLSN
         print *, 'Surface pressure, density:', P(NMP), rho(NMP)
         print *, 'Effective temperature:    ', T(NMP)
         print *, 'Central temperature:      ', T(1)
         !print *, L(NMP)/(T(NMP)**4*r(NMP)**2) / (CLSN/CRSN**2)
         print *, 'Final residuals:          ', PBC, MBC
         print *, ''
      end if

      return

!     Check self-consistency of solution: we assumed full convection
      IJ = 0
      do i=1, NMP
         GRADR = 3.0d0*KAPPA(i)*P(i)*L(i) / (16.0*CPI*CL*CG*CA*T(i)**4*(M(i) + 1.0d-16))
         str = "  (stable)"
         if (GRADAp(i) > GRADR) IJ = IJ+1
      enddo
      if (verbose) then
         print '(I4, " of ", I4, " meshpoints should be radiative")', IJ, NMP
      end if

      open (unit=66, file="polytrope.out")
      write (66, '("#",1P, 1X, A4, 9A12)') 'MP', 'm', 'r', 'rho', 'P', 'T', 'L', 'Q', 'QQ', 'conv.'
      do i=1, NMP
         GRADR = 3.0d0*KAPPA(i)*P(i)*L(i) / (16.0*CPI*CL*CG*CA*T(i)**4*(M(i) + 1.0d-16))
         str = "  (stable)"
         if (GRADAp(i) < GRADR) str = "  (unstable)"
         write (66, '(1P,1X,i5, 8E12.4, A12)') i, m(i), r(i), rho(i), p(i), T(i), L(i), Q(i), QQ(i), str
      end do

!      print *, (Q(NMP)-Q(1)), QQ(NMP)-QQ(1)

   end subroutine get_polytrope



!     calculate polytrope parameters: n, a, k, rho_c
   subroutine compute_polytrope_constants(mass, radius)
   use constants
   use poly32
   implicit none
   real(double), intent(in) :: mass, radius

   n = 3./2.            ! Polytropic index
   avg_rho = mass/(4./3.*CPI*radius**3)
   rho_c = avg_rho / (3.0/xi32(n_poly_mesh)**3 * dthetadxi32(n_poly_mesh))
   a = xi32(n_poly_mesh) / radius
   k = 4.0*CPI*CG/((n+1.0)*a**2) * rho_c**((n-1.0)/n) * 1.0e22
   end subroutine compute_polytrope_constants



!     Based on the values for polytrope, construct an Eggleton-type mesh point
!     distribution. This is relatively straightforward because the terms in
!     normal mesh spacing function can all be expressed in terms of known
!     values and prameters for the polytrope.
   subroutine setup_mesh_spacing(star_radius, poly_radius, xa)
   use constants
   use settings
   use poly32
   implicit none
   real(double), intent(in) :: star_radius, poly_radius, xa(9)
   real(double) :: xi, theta, dthetadxi, dq, dqi,QP(n_poly_mesh)
   real(double) :: VP, VT, VM, VR, VMF, XMC, gamma, T0, P_c
   real(double) :: logf, logT
   integer :: i, nn, nh, nl

   do i = 1, NMP
      r(i) = float(i-1)/float(NMP-1)*star_radius
   end do

   P_c = K * rho_c**n
   call prtoft(log(P_c), log(rho_c), logf, logT, xa)
   T0 = exp(logT)
   do i = 1, n_poly_mesh
      !r(i) = float(i-1)/float(NMP-1)*radius
      !xi = a*r(i)
      xi = xi32(i) * star_radius / poly_radius
      call interpolate_poly32(xi, theta, dthetadxi)

      gamma = 1.0 + 1.0/n
      !wmu = 2.0 / (1.0 + 3.0*XA(1) + 0.5*XA(2))
      !T0 = WMU * AMU*k * rho_c**(1.0/n) / BOLTZM

      VP = CT(4)*DLOG(k*rho_c**gamma*theta**(n+1.0)) + CT(5)*DLOG(theta**(n+1.0) + CT(9)/(k*rho_c**gamma))
      VP = CT(4)*DLOG(theta**(n+1.0)) + CT(5)*DLOG(theta**(n+1.0) + CT(9)/(k*rho_c**gamma))
      VT = CT(7)*DLOG(T0*theta/(T0*theta + CT(10)))
! VR and VM must go to zero at centre like r**2, m**(2/3)
      XMC = CT(6)*(1./3.*CT(8)**1.5)**(2./3.)
      VMF = XMC + (dthetadxi * CPI/a**3)**(2./3.)
      VM = DLOG(DABS(XMC/VMF))
      VR = -CT(3)*DLOG(xi**2/(a**2*CT(8)) + 1.0D0)
! Q(i) is the quantity that the meshpoints should be at equal intervals of
      QP(i) = VP + VT + VM + VR
   enddo
! Construct new mesh spacing function
   QC = (QP(n_poly_mesh) - QP(1)) / float(NMP-1)
   QQ(1) = QP(1)
   do i = 2, NMP
      QQ(i) = QQ(i-1) + QC
   end do
! Now find xi as a function of Q
   do i=1, NMP
      ! Binary search
      nh = n_poly_mesh
      nl = 1
      do while (nh > nl+1)
         nn = nl + (nh - nl) / 2
         if (QP(nn) < QQ(i)) then      ! n is new upper bound
            nh = nn
         else                          ! n is new lower bound
            nl = nn
         end if
      end do 
      !do nl = 1, n_poly_mesh-1
      !   if (QP(nl+1) < QPP(i)) exit;
      !end do
      !print *, QPP(i), QP(nl), QP(nh) 

      dq = min(1.0d0, (QQ(i) - QP(nl)) / (QP(nh) - QP(nl)))
      dqi = 1.0d0 - dq
      xi = dqi*xi32(nl) + dq*xi32(nh)

      r(i) = xi/a * star_radius/poly_radius
   end do
   end subroutine setup_mesh_spacing



   subroutine integrate_polytrope(mass, xa)
   use constants
   use settings
   use poly32
   use eostate_types
   use equation_of_state
   implicit none
   real(double), intent(in) :: mass, xa(9)
   real(double) :: xi, theta, dthetadxi, Ubind
   integer :: i
   real(double) :: dm(NMP), CPp(NMP), DELTAp(NMP)
   real(double) :: Enuc(NMP), Egrav(NMP), Lgrav, Lscale
   real(double) :: VP, VT, VM, VR, VMF, XMC
   real(double) :: logT, logf
   type(eostate) :: eos
   type(abundance) :: abund

   m(0) = 0.0d0
   r(0) = 0.0d0
   L(0) = 0.0d0
   Inert(0) = 0.0d0
   Ubind= 0.0d0
   Lnuc = 0.0d0
   do i = 1, NMP
      xi = a*r(i)
      call interpolate_poly32(xi, theta, dthetadxi)
      rho(i) = rho_c * theta**n
      P(i)   = k * rho(i)**(1.0+1.0/n)
      dm(i)  = rho(i) * 4./3.*CPI * (r(i)**3 - r(i-1)**3)
      m(i)   = m(i-1) + dm(i)
      Inert(i) = Inert(i-1) + 2./3. * r(i)**2*dm(i)
      if (i>1) Ubind = Ubind + 4.0d55*CPI*CG * m(i)*dm(i) / r(i)

      !call prtoft(log(P(i)), log(rho(i)), logf, logT)
      !rho(i) = rho(i) * (ne + nio) / (neo + nio)
      call prtoft(log(P(i)), log(rho(i)), logf, logT, xa)
      call statef(logf,logT, xa(:), abund, eos)
      call nucrat (logt, abund, eos )

      !mu = 2.0 / (1.0 + 3.0*XA(1) + 0.5*XA(2))
      !print *, i, TT / ( P(i) / (CR * rho(i) * mu) )
      T(i)      = eos%T
      lnf(i)    = logf
      CPp(i)    = eos%SCP
      GRADAp(i) = eos%GRADA
      DELTAp(i) = eos%DELTA
      kappa(i)  = eos%FK
      Enuc(i)   = eos%EX
      Lnuc      = Lnuc + Enuc(i)*dm(i)
! Mesh spacing equation
      VP = CT(4)*LOG(P(i)) + CT(5)*LOG(P(i) + CT(9))
      VT = CT(7)*LOG(T(i)/(T(i) + CT(10)))
! VR and VM must go to zero at centre like r**2, m**(2/3)
      XMC = CT(6)*(4./3.*rho_c*CT(8)**1.5)**(2./3.)
      VMF = XMC + M(i)**(2./3.)
      VM = LOG(ABS(XMC/VMF))
      VR = -CT(3)*LOG(R(i)**2/CT(8) + 1.0D0)
! Q(i) is the quantity that the meshpoints should be at equal intervals of
      Q(i) = VP + VT + VM + VR
   end do

   !print *, Lnuc, L(NMP)

!  Integrate luminosity equation
!  We assume that the energy is being produced by homologous contraction
!  This gives eps_grav = -3/5*cp*T * Rdot/R.
!  Integrating this expression fixes L up to the overall scale factor
!  Rdot/R, which can then be determined from the known luminosity at the
!  surface.
   L(1) = 0.0
   L(NMP) = 4.0*CPI*R(NMP)**2 * 0.25*CA*CL*T(NMP)**4 * 1.0d-11
   Lgrav = 0.0
   do i = 2, NMP
      Egrav(i) = 0.6 * CPp(i) * T(i)
      Lgrav = Lgrav + Egrav(i) * dm(i)
   end do
   Lscale = (L(NMP) - Lnuc) / Lgrav
   Egrav = Egrav * Lscale
   do i = 2, NMP-1
      L(i) = L(i-1) + (Egrav(i) + Enuc(i)) * dm(i)
   end do

!  Kelvin-Helmholtz timescale
   TKH = 0.5d-33 * Ubind / L(NMP)

!  Boundary conditions
   PBC = (eos%PG + 0.5*eos%PR) - 2./3. * 1.0d11*CG*m(NMP)/(r(NMP)**2*eos%FK)
   MBC = mass - m(NMP)

   end subroutine integrate_polytrope



   subroutine generate_polytrope_model(mass_in_msun, Tc, xa, h, dh, hnuc)
   use constants
   use settings
   use atomic_data
   use poly32
   use opacity
   use mesh, only:nm, nvar
   use indices
   implicit none
   real(double), intent(in) :: mass_in_msun, Tc, xa(9)
   real(double), intent(out) :: h(nvar,nm), dh(nvar,nm), hnuc(50, nm)
   real(double) :: mass, radius, rr, dm, dr, dtdr
   integer :: i, kk

   ! Read in polytrope structure
   call read_poly32

   ! Initialise opacity table
   KOP = 1
   KION = 5
   call load_opacity(kop, czs)

   ! Set coefficients for mesh spacing function
   CT = 0.0d0
   !CT(1:10) = (/ 0.0E+00, 0.0E+00, 5.0E-02, 5.0E-02, 0.15E+0, 2.0E-02, 0.45E+0, 1.0E-04, 1.0E+15, 2.0E+04 /)
   CT(3) = 5.0d-02
   CT(4) = 5.0d-02
   CT(5) = 0.15d0
   CT(6) = 2.0d-02
   CT(7) = 0.45d0
   CT(8) = 1.0d-4
   CT(9) = 1.0d15
   CT(10) = 2.0d4

   ! Convert mass to code units
   mass = mass_in_msun * CMSN

   ! Crude initial guess for the radius
   radius = 8.0d0 * CRSN * mass / CMSN
   radius = min(400.d0 * CRSN, radius)
   radius = max(3.d0 * CRSN, radius)

   ! Initial guess for corrections
   dm = 0.d0!mass * 7.0d-2
   dr = radius * 1.0d-2
   !radius = radius - dr

   print *, 'Generating polytropic starting model for mass ', mass_in_msun
   print *, 'Initial guess: ', (mass+dm)/CMSN, (radius+dr)/CRSN
   do i=1, 20
      call get_polytrope(mass, radius, xa, dm, dr)
      dtdr = T(1)
      call get_polytrope(mass, radius + 1.0d-6, xa, dm, dr)
      dtdr = (T(1) - dtdr) / 1.0d-6;
      rr = (Tc - T(1)) / dtdr
      if (radius + rr < 0.0d0) rr = -radius * 1.0d-1
      if (radius + rr > 1.1*radius) rr = radius * 1.0d-1
      radius = radius + rr
      call get_polytrope(mass, radius, xa, dm, dr)

      print *, 'Polytrope mass and radius:', (mass+dm)/CMSN, (radius+dr)/CRSN
      print *, 'Stellar mass and radius:  ', m(NMP)/CMSN, r(NMP)/CRSN
      print *, 'Luminosity (total, nuc.): ', L(NMP)/CLSN, Lnuc/CLSN
      print *, 'Central pressure, density:', P(1), rho(1)
      print *, 'Surface pressure, density:', P(NMP), rho(NMP)
      print *, 'Effective temperature:    ', T(NMP)
      print *, 'Central temperature:      ', T(1)

      if (abs((T(1) - Tc)/Tc) < 1.0d-6) exit
   end do

   ! Convert output to format used by the evolution code
   H(:,:) = 0.0d0
   do i = 1, NMP
      kk = NMP + 1 - i
      H(VAR_LNF,kk)   = lnf(i)
      H(VAR_LNT,kk)   = log(T(i))
      H(VAR_MASS,kk)  = m(i)
      H(VAR_QK,kk)    = QC
      H(VAR_LNR,kk)   = 0.5*log(r(i)**2 + CT(8))
      H(VAR_LUM,kk)   = L(i)
      H(VAR_INERT,kk) = Inert(i)
      H(VAR_OMEGA,kk) = 1.0d6

      H(VAR_H1, kk)   = XA(1)
      H(VAR_HE4, kk)  = XA(2)
      H(VAR_C12, kk)  = XA(3)
      H(VAR_N14, kk)  = XA(4)
      H(VAR_O16, kk)  = XA(5)
      H(VAR_NE20, kk) = XA(6)
      H(VAR_MG24, kk) = XA(7)
      H(VAR_SI28, kk) = XA(8)
      H(VAR_FE56, kk) = XA(9)
   end do

   dh = 0.0d0
   hnuc = 0.0d0

   end subroutine generate_polytrope_model



   subroutine generate_starting_model(mass_in_msun, h, dh, hnuc, sm,dty,age,per,bms,ecc,p1,enc,kh,kp,jmod,jb,jn,jf)
   use constants
   use settings
   use atomic_data
   use poly32
   use opacity
   use mesh, only:nm, nvar
   use indices
   implicit none
   real(double), intent(in) :: mass_in_msun
   integer, intent(out) :: kh,jn,jf, jmod,kp,jb
   real(double), intent(out) :: sm,dty,age,per,bms,ecc,p1,enc
   real(double), intent(out) :: h(nvar,nm), dh(nvar,nm), hnuc(50, nm)
   real(double) :: Tc
   real(double) :: ch0, che, vma, xa(9)

   ! Set coefficients for mesh spacing function
   CT = 0.0d0
   !CT(1:10) = (/ 0.0E+00, 0.0E+00, 5.0E-02, 5.0E-02, 0.15E+0, 2.0E-02, 0.45E+0, 1.0E-04, 1.0E+15, 2.0E+04 /)
   CT(3) = 5.0d-02
   CT(4) = 5.0d-02
   CT(5) = 0.15d0
   CT(6) = 2.0d-02
   CT(7) = 0.45d0
   CT(8) = 1.0d-4
   CT(9) = 1.0d15
   CT(10) = 2.0d4

   ! We want a central temperature of ~1e6K, where all nuclear reactions switch
   ! off.
   Tc = 1.0d6

   ! ZAMS composition, from init.dat
   ch0 = ch
   if (czs < 0.0d0) czs = 0.02d0
   if (ch<0.0d0) ch0 = 0.76d0 - 3.0d0*czs
   che = 1.0d0 - ch0 - czs
   xa(1) = ch0*cbn(1)/can(1)
   xa(2) = che*cbn(2)/can(2)
   xa(3) = cc*czs*cbn(3)/can(3)
   xa(4) = cn*czs*cbn(4)/can(4)
   xa(5) = co*czs*cbn(5)/can(5)
   xa(6) = cne*czs*cbn(6)/can(6)
   xa(7) = cmg*czs*cbn(7)/can(7)
   xa(8) = csi*max(czs, 1.0d-4)*cbn(8)/can(8)
   xa(9) = cfe*max(czs, 1.0d-4)*cbn(9)/can(9)
   vma = sum(xa(1:9))
   xa(1:9) = xa(1:9) / vma

   call generate_polytrope_model(mass_in_msun, Tc, xa, H, dH, Hnuc)

   sm   = m(NMP) / CMSN
   dty  = 0.1*TKH/CSY
   age  = 0.0d0
   per  = 1.0d6
   bms  = 0.1 + sm
   ecc  = 0.0d0
   p1   = 1.0d6
   enc  = 0.0d0
   kh   = NMP
   kp   = 2000
   jmod = 0
   jb   = 1
   jn   = 24
   jf   = 0

   end subroutine generate_starting_model
end module polytrope
