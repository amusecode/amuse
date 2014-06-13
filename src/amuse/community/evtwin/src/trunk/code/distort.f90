! Functions for calculating the shape of a distorted star, either due to tides or rotation.
! Based on a multipole expansion of the gravitational potential of the distorted star, as described
! by Kopal (1959) and more recently Landin & al. (2009), A&A 494.
!
! Main differences:
!  1. Do not neglect the contribution of the surface quadrupole momement in the interior.
!     It can become important deep in the stellar interior.
!  2. Due to (1), there are deviations from Radau's equations as given in the literature.
!     In particular, the equation for the l=2 mode is different for rotation and tides.
!  3. Do not assume that rotation is uniform, but assume shellular rotation.
!     This means that the density is also not constant on isobars (TODO).
!
! Simplifying assumptions:
!  1. The tidal potential of the companion is that of a point mass. It should be doable (but cumbersome) to do
!     the full treatment, but it may not be worth it. Certainly the next order (using components of the inertia
!     tensor) should be tried first.
!  2. The expansion is linear in the coefficients of spherical harmonics; higher order terms are ignored.
!  3. The rotation axis is perpendicular to the orbital plane (rotational and orbital angular momentum is aligned).
!     The coordinate system corotates with the orbit, coordinate axes are chosen such that the rotation axis is along
!     the z-direction and the line connecting the two stars is along the x-axis.
!  4. Stars are made up of perfect fluid, so there is no tidal lag and no m != 0 modes.
!     The multipole expansion then reduces to an expansion in Legendre polynomials.
module distortion
   use real_kind
   use constants, only: CPI

   ! Instead of integrating over the full solid angle, we can get away with integrating over only one quarter, since
   ! the rotational perturbation has xy-mirror symmetry and the tidal perturbation has xz-mirror symmetry.
   ! Note that rotation alone would also have xz and yz symmetries, but the yz symmetry is broken by the tides.
   real(double), parameter, private :: theta_min       = 0.0
   real(double), parameter, private :: theta_max       = CPI / 2.0
   real(double), parameter, private :: phi_min         = 0.0
   real(double), parameter, private :: phi_max         = CPI
   real(double), parameter, private :: symmetry_factor = 4.0

   ! Number of steps to take in the phi, theta directions
   ! NB: in the theta direction, we should perhaps use steps that are equidistant in cos(theta) rather than theta
   ! itself...
   integer, parameter, private :: n_theta     = 101
   integer, parameter, private :: n_phi       = 101
   real(double), parameter, private :: dtheta = (theta_max - theta_min) / n_theta
   real(double), parameter, private :: dphi   = (phi_max - phi_min) / n_phi

   ! Values of theta, phi in each direction
   real(double), private, save :: vtheta(n_theta)
   real(double), private, save :: vphi(n_phi)

   ! TODO: test whether it is faster to pre-calculate sine and cosines...

contains

   subroutine initialise_distortion_tables
   implicit none
   integer :: n

   forall(n=1:n_theta) vtheta(n) = (n-1) * (theta_max - theta_min) / (n_theta - 1)
   forall(n=1:n_phi) vphi(n) = (n-1) * (phi_max - phi_min) / (n_phi - 1)

   do n=1, n_theta
      !print *,  n, vtheta(n), asin( real(n-1) / (n_theta-1) )
   end do

   end subroutine initialise_distortion_tables

! ------------------------------------------------------------------------------
! legendrep
! Calculate a Legendre polynomial of degree l, for parameter x.
! ------------------------------------------------------------------------------
! Input:
!  l     - degree of the polynomial (we only implement up to l=5)
!  x     - parameter value
! Returns:
!  The value of the Legendre polynomial
! ------------------------------------------------------------------------------
   elemental function legendrep(l, x)
   implicit none
   integer, intent(in) :: l
   real(double), intent(in) :: x
   real(double) :: legendrep

   legendrep = 1.
   select case (l)
      case (1)
         legendrep = x
      case (2)
         legendrep = 0.5*(3.*x**2 - 1.)
      case (3)
         legendrep = 0.5*(5.*x**3 - 3.*x)
      case (4)
         legendrep = 0.125*(34.*x**4 - 30.*x**2 + 3.)
      case (5)
         legendrep = 0.125*(63.*x**5 - 70.*x**3 + 15.*x)
   end select
   end function legendrep
   


! ------------------------------------------------------------------------------
! dlegendrep
! Calculate the derivative of a Legendre polynomial of degree l, for parameter x.
! ------------------------------------------------------------------------------
! Input:
!  l     - degree of the polynomial (we only implement up to l=4)
!  x     - parameter value
! Returns:
!  The value of the derivative of the Legendre polynomial
! ------------------------------------------------------------------------------
   elemental function dlegendrep(l, x)
   implicit none
   integer, intent(in) :: l
   real(double), intent(in) :: x
   real(double) :: dlegendrep

   dlegendrep = 0.
   select case (l)
      case (1)
         dlegendrep = 1.
      case (2)
         dlegendrep = 3.*x
      case (3)
         dlegendrep = 1.5*(5.*x*x - 1.)
      case (4)
         dlegendrep = 2.5*(7.*x*x*x - 3.*x)
   end select
   end function dlegendrep



! ------------------------------------------------------------------------------
! intradau
! Integrate Radau's equation for the logarithmic derivative of the spherical
! harmonic.
! ------------------------------------------------------------------------------
! Input:
!  nm          - the number of mesh points, ordered surface (1) to centre (nm)
!  l           - the order of the spherical harmonic
!  rhoratio(:) - ratio of density to average interior density
!  radius0(:)  - average radius of each shell [1e11 cm]
!  dradius0(:) - average thickness of each shell [1e11 cm]
! Ouput:
!  eta(:)      - logarithmic derivative of spherical harmonic at each shell
! ------------------------------------------------------------------------------
   pure subroutine intradau(nm, l, rhoratio, radius0, eta)
   implicit none
   integer, intent(in)        :: nm, l
   real(double), intent(in)   :: rhoratio(nm), radius0(nm)
   real(double), intent(out)  :: eta(nm)
   real(double)               :: dr0
   integer :: n

   ! Central boundary condition value
   eta(nm) = l - 2

   ! Integrate from centre to surface
   ! FIXME: use something better than plain forward Euler, say Runge-Kutta...
   !        The main difficulty is that we'll also need to interpolate the density for that.
   do n = nm-1, 1, -1
      dr0 = radius0(n) - radius0(n+1)
      eta(n) = eta(n+1) + deta(eta(n+1), radius0(n+1), rhoratio(n+1), l) * dr0
   end do

   contains

      ! Right hand side of Radau's equation, for the integration function
      elemental function deta(eta, r0, rhoratio, l)
      implicit none
      real(double), intent(in) :: eta, r0, rhoratio
      integer, intent(in) :: l
      real(double) :: deta
      deta = eta*(eta-1.) + 6.*rhoratio*(1.+eta) - l*(l+1)
      deta = -deta / r0
      end function deta

   end subroutine intradau



! ------------------------------------------------------------------------------
! compute_envelope_integral
! Compute the envelope integral terms that enter in the expression for the
! spherical harmonics.
! ------------------------------------------------------------------------------
! Input:
!  nm          - the number of mesh points, ordered surface (1) to centre (nm)
!  l           - the order of the spherical harmonic
!  rho(:)      - the avarega density on an isobar [g/cm^3]
!  radius0(:)  - average radius of each shell [1e11 cm]
!  alpha(:)    - coefficient of spherical harmonic
!  eta(:)      - logarithmic derivative of spherical harmonic at each shell
! Output:
!  sigma(:)    - envelope integral
! ------------------------------------------------------------------------------
   pure subroutine compute_envelope_integral(nm, l, rho, radius0, alpha, eta, sigma)
   use constants, only: CG, CPI4
   implicit none
   integer, intent(in)       :: nm, l
   real(double), intent(in)  :: rho(nm), radius0(nm), alpha(nm), eta(nm)
   real(double), intent(out) :: sigma(nm)
   real(double) :: dr
   integer :: n

   dr = radius0(2) - radius0(1)
   sigma(1) = rho(1) * radius0(1)**(1 - l)*alpha(1)*(2.-l+eta(1)) * dr

   ! Integrate surface->centre (or gridpoint to surface)
   do n=2, nm
      dr = radius0(n-1) - radius0(n)
      sigma(n) = sigma(n-1) + rho(n) * radius0(n)**(1 - l)*alpha(n)*(2.-l+eta(n)) * dr
   end do

   sigma = sigma * CG * CPI4 / (2.*l + 1)
   end subroutine compute_envelope_integral



! ------------------------------------------------------------------------------
! compute_interior_integral
! Compute the interior integral terms that enter in the expression for the
! spherical harmonics.
! Note that l=0 gives the interior mass
! ------------------------------------------------------------------------------
! Input:
!  nm          - the number of mesh points, ordered surface (1) to centre (nm)
!  l           - the order of the spherical harmonic
!  rho(:)      - the avarega density on an isobar [g/cm^3]
!  radius0(:)  - average radius of each shell [1e11 cm]
!  alpha(:)    - coefficient of spherical harmonic
!  eta(:)      - logarithmic derivative of spherical harmonic at each shell
! Output:
!  mu(:)       - envelope integral
! ------------------------------------------------------------------------------
   pure subroutine compute_interior_integral(nm, l, rho, radius0, alpha, eta, mu)
   use constants, only: CPI4
   implicit none
   integer, intent(in)       :: nm, l
   real(double), intent(in)  :: rho(nm), radius0(nm), alpha(nm), eta(nm)
   real(double), intent(out) :: mu(nm)
   real(double) :: dr
   integer :: n

   dr = radius0(nm) - radius0(nm-1)
   mu(nm) = rho(nm) * radius0(nm)**(2 + l)*alpha(nm)*(3.+l+eta(nm)) * dr

   ! Integrate centre -> surface
   do n=nm-1, 1, -1
      dr = radius0(n) - radius0(n+1)
      mu(n) = mu(n+1) + rho(n) * radius0(n)**(2 + l)*alpha(n)*(3.+l+eta(n)) * dr
   end do

   mu = mu * CPI4 / (2.*l + 1)
   end subroutine compute_interior_integral



! ------------------------------------------------------------------------------
! calculate_distortion_coefficients
! Find the coefficients for the spherical harmonics
! ------------------------------------------------------------------------------
! Input:
!  nm          - the number of mesh points, ordered surface (1) to centre (nm)
!  rho(:)      - the avarega density on an isobar [g/cm^3]
!  radius(:)   - volume radius of the isobar [1e11 cm]
!  r0_guess(:) - Initial guess for the average radius [1e11 cm], could be volume radius
!  mass(:)     - enclosed mass for isobar [1e33 g]
!  omega(:)    - rotation rate for each isobar (shellular rotation) [s^-1]
!  m2          - the mass of the secondary star [1e33 g]
!  a_orb       - the orbital separation between the two stars [1e11 cm]
! Output:
!  r0(:)       - the average radius of each isobar
!  a_rot2(:)   - coefficient of the spherical harmonic for rotation (l=2)
!  a_tide(:,:) - coefficients for tidal terms (l=2,3,4)
!  sigma_rot(:)- exterior integral for rotation (l=2)
!  mu_rot(:)   - interior integral for rotation (l=2)
!  sigma_tide(:,:)- exterior integral for tides (l=2,3,4)
!  mu_tide(:,:)- interior integral for tides (l=2,3,4)
! ------------------------------------------------------------------------------
   pure subroutine calculate_distortion_coefficients(nm, rho, r, m, omega, m2, a_orb, r0, a_rot2, a_tide,&
                                                                                      sigma_rot, mu_rot, sigma, mu)
   use constants, only: CG, CPI4, C3RD
   implicit none
   integer, intent(in) :: nm
   real(double), intent(in) :: rho(nm), r(nm), m(nm), omega(nm)
   real(double), intent(in) :: m2, a_orb
   real(double), intent(out) :: r0(nm), a_rot2(nm), a_tide(nm,2:4), mu_rot(nm), mu(nm,2:4), sigma_rot(nm), sigma(nm,2:4)
   real(double) :: q, rho_ratio(nm), eta(nm, 2:4)
   real(double) :: epsilon, delta, a, b, c, d, a2, b2, c2, d2, f3, f
   integer :: iter, n, j

   ! Initialise with default values: undistorted star
   r0 = r
   a_rot2 = 0
   a_tide = 0
   sigma_rot = 0.0
   sigma = 0.0

   ! We are looking for r0, using r as an initial guess (FIXME: we can do better by taking r0 from the previous iteration).
   ! In general we will need to iterate a number of times for the solution to converge.
   do iter=1, 10
      ! Calculate coefficients for Radau's equation
      rho_ratio(:) = rho(:) / ( m(:) / (CPI4*C3RD * r0(:)**3) )
      rho_ratio(nm) = 1.0

      ! Integrate Radau's equations
      ! For now, there's no difference between the l=2 modes for rotation and for tides
      call intradau(nm, 2, rho_ratio, r0, eta(:,2))
      call intradau(nm, 3, rho_ratio, r0, eta(:,3))
      call intradau(nm, 4, rho_ratio, r0, eta(:,4))

      ! Calculate envelope integrals.
      ! We only do this at the second iteration, since values are undefined for the first iteration
      if (iter > 1) then
         call compute_envelope_integral(nm, 2, rho, r0, a_rot2,      eta(:,2), sigma_rot)
         if (m2 > 0.0) then
            call compute_envelope_integral(nm, 2, rho, r0, a_tide(:,2), eta(:,2), sigma(:,2))
            call compute_envelope_integral(nm, 3, rho, r0, a_tide(:,3), eta(:,3), sigma(:,3))
            call compute_envelope_integral(nm, 4, rho, r0, a_tide(:,4), eta(:,4), sigma(:,4))
         end if
      end if

      delta = 0.0

      do n=1, nm
         ! Calculate expansion coefficients

         ! Rotation
         a_rot2(n) = -5./(2. + eta(n, 2)) * (omega(n)**2/3. + sigma_rot(n)) * r0(n)**3/(CG*m(n));

         ! Tides
         if (m2 > 0.0) then
            q = m2/m(n)
            forall (j=2:4) &
               a_tide(n,j) = (2.*j+1.)/(j + eta(n,j)) * q * (r0(n)/a_orb)**(j+1) * (1.0 + sigma(n,j)*a_orb**(1+j)/(CG*m2))
         end if

         ! Calculate average radius from volume radius and the shape of the potential. Test for convergence
         a = -a_rot2(n)
         b =  a_tide(n, 2)
         c =  a_tide(n, 3)
         d =  a_tide(n, 4)
         a2 = a*a;
         b2 = b*b;
         c2 = c*c;
         d2 = d*d;
         f3 = 1.0 + 3./5.*(a2+a*b+b2) + 3.*c2/7.                                                   &
            - 1./35.*(2*a*a2 + 3*b*a2 - 3.*a*b2-2.*a*c2-2*b*b2 - 6.*a*b*c - 6.*b2*d - 4.*b*c2)     &
            + 10./231.*(a*d2 + 2*b*d2)                                                             &
            + 9./140.*a2*d + 6.*c2*d/77 + 18./1001.*d*d2 + d2/3.
         f = f3**C3RD

         ! Change to radius in this iteration
         epsilon = abs(r(n)/f - r0(n))
         if (epsilon > delta) delta = epsilon

         ! Update average radius for next iteration
         r0(n) = r(n) / f
      end do
      !print *, iter, delta
      if (delta < 1.0e-12) exit
   end do

   ! Calculate interior integrals
   mu_rot = 0.0
   mu = 0.0
   !call compute_interior_integral(nm, 2, rho, r0, a_rot2,      eta(:,2), mu_rot)
   !call compute_interior_integral(nm, 2, rho, r0, a_tide(:,2), eta(:,2), mu(:,2))
   !call compute_interior_integral(nm, 3, rho, r0, a_tide(:,3), eta(:,3), mu(:,3))
   !call compute_interior_integral(nm, 4, rho, r0, a_tide(:,4), eta(:,4), mu(:,4))
   end subroutine calculate_distortion_coefficients



! ------------------------------------------------------------------------------
! calculate_effective_gravity
! Calculate the average effective gravity <g> and its inverse <1/g> over an
! isobar, as well as the surface area of the isobar.
! To rotational component is invariant under reflection in the xy, xz and yz
! planes, but the tidal component is asymmetric along the x-axis. This means
! we can break up the integration domain in 4 equal parts and integrate the
! azimuthal angle phi from 0 to pi and the polar angle theta from 0 to pi/2.
! ------------------------------------------------------------------------------
! Input:
!  nm            - the number of mesh points, ordered surface (1) to centre (nm)
!  radius(:)     - volume radius of the isobar [1e11 cm]
!  mass(:)       - enclosed mass for isobar [1e33 g]
!  omega(:)      - rotation rate for each isobar (shellular rotation) [s^-1]
!  m2            - the mass of the secondary star [1e33 g]
!  a_orb         - the orbital separation between the two stars [1e11 cm]
!  r0(:)         - the average radius of each isobar
!  a_rot2(:)     - coefficient of the spherical harmonic for rotation (l=2)
!  a_tide(:,:)   - coefficients for tidal terms (l=2,3,4)
!  sigma_rot2(:) - envelope integral for rotation (l=2)
!  gamma_rot2(:) - interior integral for rotation (l=2)
!  sigma_tide(:) - envelope integral for tides (l=2,3,4)
!  gamma_tide(:) - interior integral for tides (l=2,3,4)
! Output:
!  avgg(:)       - effective gravity over an isobar (cm / s^2)
!  invg(:)       - inverse effective gravity over an isobar (cm / s^2)
!  spsi(:)       - surface area of an isobar (1e22 cm^2)
!  fp(:)         - correction factor fp [-]
!  ft(:)         - correction factor ft [-]
! ------------------------------------------------------------------------------
   pure subroutine calculate_effective_gravity(nm, ra, ma, omega, m2, a_orb, r0a, a_rot2, a_tide,&
                                               sigma_rot2, gamma_rot2, sigma_tidea, gamma_tidea, avgg, invg, spsi, fp, ft)
   use constants, only: CG, CPI4
   implicit none
   integer, intent(in) :: nm
   real(double), intent(in) :: ra(nm), ma(nm), omega(nm)
   real(double), intent(in) :: m2, a_orb
   real(double), intent(in) :: a_rot2(nm), a_tide(nm,2:4), r0a(nm)
   real(double), intent(in) :: sigma_rot2(nm), gamma_rot2(nm), sigma_tidea(nm,2:4), gamma_tidea(nm,2:4)
   real(double), intent(out) :: avgg(nm), invg(nm), spsi(nm), fp(nm), ft(nm)
   real(double) :: a_orb2, omega2, theta, phi, cos_chi, rr, r2, r3, r4, dspsi, dpsi_dr, dpsi_dtheta, dpsi_dphi, fac, facr, gg
   real(double) :: sigma_rot, mu_rot, sigma_tide(2:4), mu_tide(2:4), a2, a(2:4), r, m, r0
   real(double) :: cos_theta, sin_theta, dcos_theta, cos_phi, sin_phi
   integer :: n, nt, np


   a_orb2 = a_orb*a_orb
   do n = 1, nm
      spsi(n) = 0.0
      avgg(n) = 0.0
      invg(n) = 0.0
      omega2 = omega(n)**2

      ! Copy integrals and coefficients to local variables, oterwise we keep hammering the processor cache
      sigma_rot       = sigma_rot2(n)
      mu_rot          = gamma_rot2(n)
      sigma_tide(2:4) = sigma_tidea(n, 2:4)
      mu_tide(2:4)    = gamma_tidea(n, 2:4)
      a2              = a_rot2(n)
      a(2:4)          = a_tide(n, 2:4)
      r0              = r0a(n)
      m               = ma(n)
      r               = ra(n)
      do nt=1, n_theta
         theta = (nt-1)*dtheta
         cos_theta  = cos(theta)
         sin_theta  = sin(theta)
         dcos_theta = cos_theta - cos(theta+dtheta)
         do np=1, n_phi
            phi = (np-1)*dphi
            cos_phi = cos(phi)
            sin_phi = sin(phi)
            cos_chi = cos_phi*sin_theta

            ! Calculate the radial coordinate of the isobar in this direction
            rr = 1.0 + a2 * legendrep(2, cos_theta) &
                     + a(2)*legendrep(2,cos_chi) + a(3)*legendrep(3,cos_chi) + a(4)*legendrep(4,cos_chi)
            rr = r0*rr
            r2 = rr*rr
            r3 = r2*rr
            r4 = r3*rr

            ! Surface element
            dspsi = r2 * dcos_theta*dphi
            spsi(n) = spsi(n) + dspsi

            ! Gradient of pseudo potential
            dpsi_dr     = CG*m/r2 - 2.*omega2*rr/3. &
                           -(-3.*CG*mu_rot/r3 + 2.*rr*sigma_rot - 2.*omega2*rr/3.) * legendrep(2, cos_theta)
            dpsi_dtheta = (CG*mu_rot/r3 + r2*sigma_rot - omega2*r2/3.) * dlegendrep(2, cos_theta) * sin_theta
            dpsi_dphi   = 0;

            ! Unroll the following loop manually; for some reason the Intel compiler does an appalling job with this loop,
            ! making it take much longer than it normally should.
            !do j = 2, 4
            !   fac  = CG*mu_tide(j)/rr**(j+1) + rr**j*sigma_tide(j) + CG*m2/a_orb*(rr/a_orb)**j
            !   facr = -(j+1)*CG*mu_tide(j)/rr**(j+2) + j*rr**(j-1)*sigma_tide(j) + j*CG*m2/a_orb2*(rr/a_orb)**(j-1)
            !   dpsi_dr     = dpsi_dr     - facr * legendrep(j, cos_chi)
            !   dpsi_dtheta = dpsi_dtheta + fac  * dlegendrep(j,cos_chi)*cos(phi)*cos(theta)
            !   dpsi_dphi   = dpsi_dphi   + fac  * dlegendrep(j,cos_chi)*sin(theta)*sin(phi)
            !end do
            fac  = CG*mu_tide(2)/r3 + r2*sigma_tide(2) + CG*m2/a_orb*(rr/a_orb)**2
            facr = -(2+1)*CG*mu_tide(2)/r4 + 2*rr*sigma_tide(2) + 2*CG*m2/a_orb2*(rr/a_orb)**(2-1)
            dpsi_dr     = dpsi_dr     - facr * legendrep(2, cos_chi)
            dpsi_dtheta = dpsi_dtheta + fac  * dlegendrep(2,cos_chi)*cos_phi*cos_theta
            dpsi_dphi   = dpsi_dphi   + fac  * dlegendrep(2,cos_chi)*sin_theta*sin_phi

            fac  = CG*mu_tide(3)/r4 + r3*sigma_tide(3) + CG*m2/a_orb*(rr/a_orb)**3
            facr = -(3+1)*CG*mu_tide(3)/(rr*r4) + 3*r2*sigma_tide(3) + 3*CG*m2/a_orb2*(rr/a_orb)**(3-1)
            dpsi_dr     = dpsi_dr     - facr * legendrep(3, cos_chi)
            dpsi_dtheta = dpsi_dtheta + fac  * dlegendrep(3,cos_chi)*cos_phi*cos_theta
            dpsi_dphi   = dpsi_dphi   + fac  * dlegendrep(3,cos_chi)*sin_theta*sin_phi

            fac  = CG*mu_tide(4)/(rr*r4) + r4*sigma_tide(4) + CG*m2/a_orb*(rr/a_orb)**4
            facr = -(4+1)*CG*mu_tide(4)/(r2*r4) + 4*r3*sigma_tide(4) + 4*CG*m2/a_orb2*(rr/a_orb)**(4-1)
            dpsi_dr     = dpsi_dr     - facr * legendrep(4, cos_chi)
            dpsi_dtheta = dpsi_dtheta + fac  * dlegendrep(4,cos_chi)*cos_phi*cos_theta
            dpsi_dphi   = dpsi_dphi   + fac  * dlegendrep(4,cos_chi)*sin_theta*sin_phi

            ! Gravity
            ! Avoid a singularity at theta=0, at which point the phi-dependent terms vanish analytically
            if (nt > 1) then
               gg = sqrt(dpsi_dr*dpsi_dr + dpsi_dtheta*dpsi_dtheta/r2 + dpsi_dphi*dpsi_dphi/(rr*sin_theta)**2)
            else
               gg = sqrt(dpsi_dr*dpsi_dr + dpsi_dtheta*dpsi_dtheta/r2)
            end if
            avgg(n) = avgg(n) + dspsi * gg
            invg(n) = invg(n) + dspsi / gg
         end do
      end do
      avgg(n) = avgg(n) / spsi(n)
      invg(n) = invg(n) / spsi(n)
      spsi(n) = spsi(n) * symmetry_factor

      ! Calculate correction factors fp and ft
      r2 = r*r
      fp(n) = CPI4 * r2*r2 / (CG*m * spsi(n) * invg(n))
      ft(n) = (CPI4*r2 / spsi(n))**2 / (avgg(n)*invg(n))

      ! Due to rounding errors, fp or ft can be > 1 (just about), which is unphysical
      fp(n) = min(fp(n), 1.0d0)
      ft(n) = min(ft(n), 1.0d0)

      ! Convert to CGS units
      avgg(n) = avgg(n) * 1.0e11
      invg(n) = invg(n) / 1.0e11
   end do
   end subroutine calculate_effective_gravity

end module distortion
