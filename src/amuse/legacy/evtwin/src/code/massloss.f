      module massloss
      use constants

! Return codes used by the different functions
! No problems encountered
      integer, parameter :: err_mdot_ok = 0
! Parameters were out of the validity range for this recipe. Proceed with caution!
      integer, parameter :: err_mdot_extrapolated_high_T =  1
      integer, parameter :: err_mdot_extrapolated_low_T  =  2
      integer, parameter :: err_mdot_extrapolated_high_L =  4
      integer, parameter :: err_mdot_extrapolated_low_L  =  8
      integer, parameter :: err_mdot_extrapolated        = 15
      integer, parameter :: err_mdot_not_used            = 16

      integer :: mdot_errno

      double precision :: METALLICITY_SCALING_EXPONENT = 0.8 

! Optionally correct de Jager for clumping (in an approximate way), but
! only for higher temperatures.
      double precision :: de_jager_clump_factor = 1.0

! Switches for the different mass loss rates

! General, to be used if none of the others are
      double precision :: multiplier_dejager = 1.0

! RGB rates
      double precision :: multiplier_schroeder = 1.0
      double precision :: multiplier_reimers = 0.0

! AGB rates
      double precision :: multiplier_vasiliadis_wood = 0.0
      double precision :: multiplier_wachter = 1.0

! A type super giants
      double precision :: multiplier_achmad = 0.0

! Hot stars
      double precision :: multiplier_vink = 1.0

! Stars close to the Eddington limit (and WR stars)
      double precision :: multiplier_kudritzki = 0.0
      double precision :: multiplier_nl = 1.0
      contains

!----------------------------------------------------!
! Helper functions for accessing bit fields in errno !
!----------------------------------------------------!

      subroutine clear_errno
      mdot_errno = 0
      end subroutine

      subroutine set_errno(error_flag)
      implicit none
      integer :: error_flag
      mdot_errno = ior(mdot_errno, error_flag)
      end subroutine

      function test_flag(flags, bitmask)
      implicit none
      logical :: test_flag
      integer :: flags, bitmask
      test_flag = iand(flags, bitmask) == bitmask
      end function


      function test_errno(error_flag)
      implicit none
      integer :: error_flag
      logical :: test_errno
      test_errno = iand(mdot_errno, error_flag) == error_flag
      end function

! calculate_mdot:
! Compute mass loss rate from stellar parameters, using various recipes
! and deciding which of those is most appropriate. The transitions
! between the various recipes are supposed to be joined together
! smoothly.
! INPUT:
!  Teff  - effective temperature, in Kelvin
!  L     - luminosity, in solar units
!  Ledd  - Eddington luminosity, in solar units
!  M     - Mass, in solar units
!  R     - Radius, in solar units
!  X     - Surface abundance of hydrogen, by mass
!  Y     - Surface abundance of helium, by mass
!  XO    - Surface abundance fraction of oxygen, by mass
!  XC    - Surface abundance fraction of carbon, by mass
!  Z0    - Initial metallicity (proxy for the iron content)
! RETURNS:
!  The mass loss rate, in solar masses per year
      function calculate_mdot( Teff, L, Ledd, M, R, X, Y, XO, XC, Z0 )
      use constants
      implicit none
      double precision, intent(in) :: Teff, L, Ledd, M, R, X, Y, XO, XC, Z0
      double precision :: calculate_mdot

! Store results of different mass loss recipes in an array, identify
! them using symbolic constants for the array index
      integer, parameter :: DEJAGER = 1
      integer, parameter :: REIMERS = 2
      integer, parameter :: SCHROEDER = 3
      integer, parameter :: WACHTER = 4
      integer, parameter :: VASILIADIS_WOOD = 5
      integer, parameter :: VINK = 6
      integer, parameter :: KUDRITZKI = 7
      integer, parameter :: ACHMAD = 8
      integer, parameter :: NL = 9
      integer, parameter :: NUM_RECIPES = 9


! Local variables and constants
      double precision, parameter :: CMSNY = CMSN/CSY ! Solar masses per year
      double precision :: logZ
      double precision :: logL, logT

! Outcome of the different mass loss recipes
      double precision :: mdot_r(NUM_RECIPES)

! Store error return codes
      integer :: recipe_errno(NUM_RECIPES)

! Different masas loss regimes
      double precision :: mdot_agb, mdot_rgb, mdot_line, mdot_luminous
      integer :: errno_agb, errno_rgb, errno_line, errno_luminous
      double precision :: mdot_cool, mdot_hot
      integer :: errno_hot, errno_cool
      double precision :: mdot
      double precision :: xx

      mdot_r = 0.0d0
      recipe_errno = 0
      calculate_mdot = 0.0d0
      logZ = log10(Z0/CZSN)
      logL = log10(L)
      logT = log10(Teff)

! Compute contributions to the mass loss from different driving
! mechanisms

! Default prescription if no other applies:
!  * de Jager et al.
      mdot_errno = err_mdot_not_used
      if (multiplier_dejager>0.0d0)
     &   mdot_r(DEJAGER) = multiplier_dejager * calc_mdot_dejager(logL, logT)
      if (test_errno(err_mdot_extrapolated_low_L)) mdot_r(DEJAGER) = 0.0
      recipe_errno(DEJAGER) = mdot_errno

! Scale the de Jager rate with the *initial* metallicity (ie, Fe
! abundance). The other mass loss recipes either have a proper
! metallicity dependence (Vink), don't depend on metallicity (Schroeder
! & Cuntz) or depend on it in some different way (presumably dust-driven
! winds).
! Because this mass loss scaling should only apply to line driven winds,
! apply it only for Teff>8000, which is roughly where Achmad&al. find
! the border of line driven winds.
      if (Teff>9500.0) 
     &   mdot_r(DEJAGER)=mdot_r(DEJAGER)*(Z0/CZSN)**metallicity_scaling_exponent

! Wind clumping correction for the radiation driven portion of the de
! Jager rate. Turn on smoothly towards higher T
      if (Teff>8000) then
         xx = min(1.0, (Teff-8000.0)/1000.0)
         mdot_r(DEJAGER) = mdot_r(DEJAGER)*(xx*de_jager_clump_factor + (1.0-xx)*1.0)
      end if

! Mass loss mechanisms on the giant branch (driven by Alfven waves(?))
!  * Classical Reimers
!  * Eggleton version of Reimers
!  * Schroeder & Cuntz (default)
      mdot_errno = err_mdot_not_used
      if (multiplier_reimers>0.0d0)
     &   mdot_r(REIMERS) = multiplier_reimers * calc_mdot_reimers(Teff, L, M, R)
      recipe_errno(REIMERS) = mdot_errno

      mdot_errno = err_mdot_not_used
      if (multiplier_schroeder>0.0d0)
     &   mdot_r(SCHROEDER) = multiplier_schroeder * calc_mdot_schroeder(Teff, L, M, R)
      recipe_errno(SCHROEDER) = mdot_errno

! FIXME: we don't relly want to include fudges like this here. The
! proper thing to do would be to
! a. Compute all mass loss rates and store their return flags. Get all
! of them to return an "interpolation factor" that is to be used if this
! rate needs to be smoothly joined to another recipe. Return flags need
! to be a bit field because we need to know if the failure is just in
! temperature or luminosity, or both.
! b. In each regime, find the best rate. This is the rate that for which
! the parameters are "within range". If there are more valid rates, then
! take the maximum value of all valid rates. If all of them are
! extrapolated, we record that information and again take the maxiumum.
! c. If there are no valid rates for this set of parameters we
! interpolate.
! c2. Determine rates on the "hot" and "cool" side of the current set of
! parameters (M, L).
! c3. On the cool side: if all available rates fail at high and low
! luminosities as well, then interpolate the RGB rate with the AGB rate.
! c4. On the hot side: the Vink rate should not fail at low or high
! luminosities, although it is forced to 0 at smaller luminosities.
! c5. Interpolate between the high T and low T solutions, like this:
!  Mdot = xh*(1-xl)*Mdot_h + (1-xh)*xl*Mdotl
! c6. If there is a fall-back rate (de Jager):
!  return max(Mdot_fallback, Mdot)
!
! In practice, where we need this smoothing is on the hot side of
! Schroeder&Cuntz and on the cool side of Vink.
!      if (test_errno(err_mdot_extrapolated_low_T)) then
!         xx = (5000.0 - log(L) * 400.0 - Teff  ) / 2500.0
!         xx = min(1.0, max(xx, 0.0))
!         mdot_r(SCHROEDER) = xx*mdot_r(SCHROEDER) + (1.0-xx)*mdot_r(DEJAGER)
!         if (Teff>5000.0) mdot_r(SCHROEDER) = mdot_r(SCHROEDER) * (Teff-5000)**(-4)
!      end if

! Mass loss mechanisms on the AGB (dust driven)
! See papers by Schroeder for critical luminosity
!  * Classical Reimers
!  * Eggleton version of Reimers
!  * Wachter & al (default)
!  * de Jager et al
!  * van Loon
!  * Vasiliadis & Wood
      mdot_errno = err_mdot_not_used
      if (multiplier_wachter>0.0d0)
     &   mdot_r(WACHTER) = multiplier_wachter * calc_mdot_wachter(logT, logL, M)
      recipe_errno(WACHTER) = mdot_errno

      mdot_errno = err_mdot_not_used
      if (multiplier_vasiliadis_wood>0.0d0)
     &   mdot_r(VASILIADIS_WOOD) = 
     &      multiplier_vasiliadis_wood * calc_mdot_vasiliadis_wood(L, M, R)
      recipe_errno(VASILIADIS_WOOD) = mdot_errno

! Mass loss from A-type supergiants (line driven)
!  * de Jager et al
!  * Achmad (default)
      mdot_errno = err_mdot_not_used
      if (multiplier_achmad>0.0d0)
     &   mdot_r(ACHMAD) = multiplier_achmad * calc_mdot_achmad(Teff, logL, logZ)
      recipe_errno(ACHMAD) = mdot_errno

! Mass loss from hot stars (line driven)
!  * de Jager et al
!  * Vink (default)
      mdot_errno = err_mdot_not_used
      if (multiplier_vink>0.0d0)
     &   mdot_r(VINK) = multiplier_vink * calc_mdot_vink(M, L, Teff, logZ)
      recipe_errno(VINK) = mdot_errno

! Mass loss close to the Eddington limit (continuum driven, WR wind)
!  * Nugis & Lamers 2002 (default)
!  * de Jager et al.
!  * Kudritzki 2002
      mdot_errno = err_mdot_not_used
      if (multiplier_kudritzki>0.0d0)
     &   mdot_r(KUDRITZKI) = multiplier_kudritzki * calc_mdot_kudritzki(logL, Teff, logZ)
      recipe_errno(KUDRITZKI) = mdot_errno

      mdot_errno = err_mdot_not_used
      if (multiplier_nl>0.0d0)
     &   mdot_r(NL) = multiplier_nl * calc_mdot_wr_nl(logL, logT, X, Y, XO, XC, logZ)
      recipe_errno(NL) = mdot_errno

! Determine mass loss rate
! FIXME: this currently fails if both rates are valid (have no errors)
      mdot_agb = max(mdot_r(WACHTER), mdot_r(VASILIADIS_WOOD))
      errno_agb = recipe_errno(WACHTER)

      mdot_rgb = max(mdot_r(REIMERS), mdot_r(SCHROEDER))
      errno_rgb = recipe_errno(SCHROEDER)
      if (errno_rgb == err_mdot_not_used) errno_rgb = mdot_r(REIMERS)

      mdot_line = max(mdot_r(VINK), mdot_r(ACHMAD))
      errno_line = recipe_errno(VINK)

      mdot_luminous = 0.0
      errno_luminous = err_mdot_not_used
! FIXME: this commented-out bit of code blatently doesn't do what it's
! supposed to do... which is decide on which of the two rates for
! luminous stars it should use.
!      if (recipe_errno(KUDRITZKI) == err_mdot_ok) then
!         errno_luminous = err_mdot_ok
!         mdot = max(mdot_r(KUDRITZKI), mdot_luminous)
!         if (mdot_r(NL) > mdot_luminous .and. recipe_errno(NL) == err_mdot_ok) then
!            mdot = max(mdot_r(NL), mdot_luminous)
!         end if
!      elseif (recipe_errno(NL) == err_mdot_ok) then
!         mdot_luminous = mdot_r(NL)
!         mdot_luminous = max(mdot_r(NL), mdot_luminous)
!      else
!         mdot_luminous = max(mdot_r(KUDRITZKI), mdot_r(NL))
!         errno_luminous = ior(recipe_errno(KUDRITZKI), recipe_errno(NL))
!         errno_luminous = iand(errno_luminous, not(err_mdot_not_used))
!      end if
      mdot_luminous = mdot_r(NL)
      errno_luminous = recipe_errno(NL)

! Cool side
      mdot_cool = 0.0
      errno_cool = err_mdot_not_used
      if (errno_rgb == err_mdot_ok) then
         mdot_cool = mdot_rgb
         errno_cool = errno_rgb
      end if
      if (errno_agb == err_mdot_ok .and. mdot_agb>mdot_cool) then
         mdot_cool = mdot_agb
         errno_cool = errno_agb
      end if
      if (test_flag(errno_agb, err_mdot_extrapolated_low_L) .and.
     &    test_flag(errno_rgb, err_mdot_extrapolated_high_L)) then
         if (mdot_rgb>mdot_agb) then
            mdot_cool = mdot_rgb
            errno_cool = errno_rgb
         else
            mdot_cool = mdot_agb
            errno_cool = errno_agb
         end if
      end if
      if (errno_cool == err_mdot_not_used) then
         mdot_cool = mdot_rgb
         errno_cool = errno_rgb
      end if 

! Hot side
      mdot_hot = 0.0
      errno_hot = err_mdot_not_used
      if (errno_line == err_mdot_ok .or. 
     &               test_flag(errno_line, err_mdot_extrapolated_low_L)) then
         mdot_hot = mdot_line
         errno_hot = errno_line
      end if
      if (errno_luminous == err_mdot_ok .and. mdot_luminous>mdot_hot) then
         mdot_hot = mdot_luminous
         errno_hot = errno_luminous
      end if

      !if (logT<3.7) write (0, *), logT, logL, mdot_cool

      mdot = 0.0
      mdot_errno = err_mdot_not_used
      if (errno_hot == err_mdot_ok) then
         mdot = mdot_hot
         mdot_errno = errno_hot
      end if
      if (errno_cool == err_mdot_ok .and. mdot_cool>mdot) then
         mdot = mdot_cool
         mdot_errno = errno_cool
      end if
      if (test_flag(errno_hot, err_mdot_extrapolated_low_T) .and.
     &    test_flag(errno_cool, err_mdot_extrapolated_high_T) ) then
! FIXME: Add interpolation information
         mdot = 0.5*(mdot_cool+mdot_hot)
         mdot_errno = err_mdot_extrapolated
      end if
      if (mdot_errno == err_mdot_not_used) then
         mdot = mdot_cool
         mdot_errno = errno_hot
      end if

      if (mdot_errno /= err_mdot_ok) then
         mdot = max(mdot_r(DEJAGER), mdot)
      end if
      !mdot = max(mdot_agb, mdot_rgb, mdot_line, mdot_luminous)

! If the adopted rate turns out to be NL but that was extrapolated, max
! it with the de Jager rate.
!      if (test_errno(err_mdot_extrapolated_low_T) .and. mdot==mdot_r(NL))
!     &   mdot = max(mdot, mdot_r(DEJAGER))

! Fall back to the default rate in case none of those above apply
!      if (mdot < 1.0d-10) mdot = max(mdot, mdot_r(DEJAGER))

! FIXME: for hot, bright low mass stars, force a high mass loss rate.
! This is to make sure the envelope of an AGB star is stripped quickly
! once its on its way to become a white dwarf.
      !if (M<1.0 .and. logL>3.7 .and. Teff>4000.0) mdot = 1.0e-4
      calculate_mdot = mdot
      end function


c ------------------------------------------------------------------------------
c CALC_MDOT_REIMERS:
c  Calculate the mass loss rate for red giants, using the mass loss
c  prescription from Reimers (1975)
c ------------------------------------------------------------------------------
c Input parameters:
c     T: Effective temperature, in K (unused, except to check validity range)
c     L: Luminosity, in solar units
c     M: Mass, in solar units
c     R: Radius, in solar units
c Returns:
c     Mass loss rate, in solar masses per year
c ------------------------------------------------------------------------------
      FUNCTION CALC_MDOT_REIMERS ( T, L, M, R )
      use constants
      implicit none
      double precision, intent(in) :: L, M, R, T
      double precision :: CALC_MDOT_REIMERS
      double precision :: lrm, Tcrit, suppress

      call clear_errno()

! At higher effecive temperatures (above say 5000K) this mass loss rate
! no longer applies. Turn it off smoothly by modulating it with a term
! that suppresses it when T is large.
! By itself the mass loss rate would diverge.
      Tcrit = min(5000.0, 5000.0 - log10(L)*400.0)
      suppress = 1.0
      if (T>Tcrit) then
         !suppress = -1.0/(Tcrit - T)
         !call set_errno(err_mdot_extrapolated_high_T)
      end if

      lrm = L*R/M * suppress

      CALC_MDOT_REIMERS = 4.0e-13*lrm
      END FUNCTION



c ------------------------------------------------------------------------------
c CALC_CRITICAL_RGB_LUMINOSITY
c Determine minimum critical luminosity for driving a dust-driven wind, see
c Schroeder & al (1999).
c Fortunately, this diverges towards high temperatures, correctly
c stating that very hot stars don't drive a dust-driven wind.
c FIXME: should use a constraint involving the dust opacity instead, see
c Wachter & al.
c ------------------------------------------------------------------------------
c Input parameters:
c     T: Effective temperature, in Kelvin
c     logM: logarithm of the Mass, in solar units
c Returns:
c     Critical luminosity, in solar units
c ------------------------------------------------------------------------------
      FUNCTION CALC_CRITICAL_RGB_LUMINOSITY ( T, logM, M )
      use constants
      implicit none
      double precision, intent(in) :: T, logM, M
      double precision :: CALC_CRITICAL_RGB_LUMINOSITY
      double precision :: Lcrit, X

      if (T<2800.0) then
         Lcrit = log10(((T-2600.0)*5.0 + 5000.0)) + logM
      else
         x = (T-2800.0)/200.0
         Lcrit = log10( (1.0-x)*6000.0 + x*1.0d4*(M-0.05)) + logM
      end if
      CALC_CRITICAL_RGB_LUMINOSITY = Lcrit
      END FUNCTION



c ------------------------------------------------------------------------------
c CALC_MDOT_SCHROEDER:
c  Calculate the mass loss rate for red giants, using the mass loss
c  prescription from Schroeder&Cuntz (2005).
c  Mass loss is driven by Alfven waves and does NOT depend on metallicity
c ------------------------------------------------------------------------------
c Input parameters:
c     T: Effective temperature, in Kelvin
c     L: Luminosity, in solar units
c     M: Mass, in solar units
c     R: Radius, in solar units
c Returns:
c     Mass loss rate, in solar masses per year
c ------------------------------------------------------------------------------
      FUNCTION CALC_MDOT_SCHROEDER ( T, L, M, R )
      use constants
      implicit none
      double precision, intent(in) :: T, L, M, R
      double precision :: CALC_MDOT_SCHROEDER
      double precision :: g_ratio, lrm, Lcrit, Tcrit, suppress

      call clear_errno()

! At higher effecive temperatures (above say 5000K) this mass loss rate
! no longer applies. Turn it off smoothly by modulating it with a term
! that suppresses it when T is large.
! By itself the mass loss rate would diverge.
      Tcrit = min(5000.0, 5000.0 - log10(L)*400.0)
      !Tcrit = 5000.0 - log10(L)*400.0
      suppress = 1.0
      if (T>Tcrit) then
         suppress = -100.0/(Tcrit - T)
         call set_errno(err_mdot_extrapolated_high_T)
      end if

      Lcrit = CALC_CRITICAL_RGB_LUMINOSITY(T, log10(M), M)
      if (log10(L) > 1.2*Lcrit) then
         call set_errno(err_mdot_extrapolated_high_L)
      end if

! Surface gravity, in solar units (ignores centrifugal term!)
      g_ratio = R**2/M
      lrm = L*R/M * suppress

      CALC_MDOT_SCHROEDER = 8.0d-14*lrm*(T/4.0D3)**3.5 * (1.0D0+g_ratio/(4.3D3))
      END FUNCTION



c ------------------------------------------------------------------------------
c CALC_MDOT_WACHTER:
c  Calculate the mass loss rate for AGB stars, from Wachter & al. (2002)
c  Mass loss is a dust-driven superwind
c ------------------------------------------------------------------------------
c Input parameters:
c     logT: log10 effective temperature, in Kelvin
c     logL: log10 luminosity, in solar units
c     M: Mass, in solar units
c Returns:
c     Mass loss rate, in solar masses per year
c ------------------------------------------------------------------------------
      FUNCTION CALC_MDOT_WACHTER ( logT, logL, M )
      use constants
      implicit none
      double precision, intent(in) :: logT, logL, M
      double precision :: CALC_MDOT_WACHTER
      double precision :: log_mdot, Lcrit, T, logM, suppress

      call clear_errno()

! Determine minimum critical luminosity for driving a dust-driven wind, see
! Schroeder & al (1999).
! Fortunately, this diverges towards high temperatures, correctly
! stating that very hot stars don't drive a dust-driven wind.
      suppress = 0.0d0
      T = 10**logT
      logM = log10(M)
      Lcrit = CALC_CRITICAL_RGB_LUMINOSITY(T, logM, M)
      if (logL < Lcrit) then
         call set_errno(err_mdot_extrapolated_low_L)
         suppress = -5.0*(logL-Lcrit)
      end if
 
! Mass loss rate based on equation (2) in the paper, where the pulsation
! period has been eliminated.
      log_mdot = -4.52 -6.81*(logT-3.41) + 2.47*(logL-4) -1.95*logM -suppress
      CALC_MDOT_WACHTER = 10**log_mdot
      END FUNCTION



c ------------------------------------------------------------------------------
c CALC_MDOT_VASILIADIS_WOOD
c  Calculate the mass loss rate for AGB stars, from Vasiliadis&Wood 1993
c  Mass loss is a pulsation driven superwind
c ------------------------------------------------------------------------------
c Input parameters:
c     L: log10 luminosity, in solar units
c     M: Mass, in solar units
c     R: Radius, in solar units
c Returns:
c     Mass loss rate, in solar masses per year
c ------------------------------------------------------------------------------
      FUNCTION CALC_MDOT_VASILIADIS_WOOD ( L, M, R )
      use constants
      implicit none
      double precision, intent(in) :: L, M, R
      double precision :: CALC_MDOT_VASILIADIS_WOOD
      double precision :: PP, max_mdot, log_mdot, vw
      double precision :: Tcrit, T

      call clear_errno()

! Maximum value of the mass loss rate, from equation (1)
      max_mdot = 1.36e-9*L

! Find Mira pulsation period, equation (4)
      PP = 10**(-2.07 + 1.94*R - 0.9*M)

! Wind velocity, equation (3), in km/s
      vw = -13.5 + 0.056*PP
      if (vw < 3.0d0 .or. L < 1.0e3) then
         call set_errno(err_mdot_extrapolated_low_L)
      end if
 
! Find mass loss rate, from equation (2) (equation 5 for masses above
! 2.5)
      log_mdot = -11.4 + 0.0125*(PP-100.0*max(M-2.5,0.0))

! Suppress V&W mass loss rate far from the Hayashi track
      Tcrit = min(5000.0, 5000.0 - log10(L)*400.0)
      !Tcrit = 5000.0 - log10(L)*400.0
      T = (L/R**2)**0.25 * 5770
      if (T>Tcrit) then
         call set_errno(err_mdot_extrapolated_high_T)
      end if

! mdot is given by the minimum of eq. (1) and eq. (5)
      CALC_MDOT_VASILIADIS_WOOD = min(max_mdot, 10**log_mdot)
      END FUNCTION



c ------------------------------------------------------------------------------
c CALC_MDOT_DEJAGER:
c  Empirical mass loss rate for luminous star, de Jager et al (1988)
c  This mass loss rate should scale with the metallicity
c ------------------------------------------------------------------------------
c Input parameters:
c     TT: log10 effective temperature/K
c     LL: luminosity
c Returns:
c     Mass loss rate, in solar masses per year
c ------------------------------------------------------------------------------
      FUNCTION CALC_MDOT_DEJAGER ( LL, TT )
      IMPLICIT REAL*8 (A-H, L-Z)
      DIMENSION A(6,6)
      DATA A/ 
     :  6.34916D0,3.41678D0,-1.08683D0, 0.13095D0,0.22427D0,0.11968D0, !J=1
     : -5.04240D0,0.15629D0, 0.41952D0,-0.09825D0,0.46591D0,0.0D0,     !J=2
     : -0.83426D0,2.96244D0,-1.37272D0, 0.13025D0,2*0.0D0,             !J=3
     : -1.13925D0,0.33659D0,-1.07493D0, 3*0.0D0,                       !J=4
     : -0.12202D0,0.57576D0, 10*0.0D0/                                 !J=5,6
      CHEB(X, I) = DCOS(I*DACOS(MIN(MAX(X, -1.0D0), 1.0D0)))

      call clear_errno()

      MDOT = 0.0D0
      IF ( TT < 3.301D0 ) MDOT_ERRNO = ERR_MDOT_EXTRAPOLATED_LOW_T
      IF ( TT > 4.799D0 ) MDOT_ERRNO = ERR_MDOT_EXTRAPOLATED_HIGH_T
      IF ( LL < 2.501D0 ) MDOT_ERRNO = ERR_MDOT_EXTRAPOLATED_LOW_L
      IF ( LL > 6.699D0 ) MDOT_ERRNO = ERR_MDOT_EXTRAPOLATED_HIGH_L
      DO I = 0, 5
         DO J = 0, 5 - I
            MDOT = MDOT + A(I + 1, J + 1)*CHEB((TT - 4.05D0)/0.75D0, I)
     :                                   *CHEB((LL - 4.60D0)/2.10D0, J)
         END DO
      END DO
      CALC_MDOT_DEJAGER = 10**(-MDOT)
      END FUNCTION



c ------------------------------------------------------------------------------
c CALC_MDOT_ACHMAD:
c Calculate the mass loss rate of Achmad, Lamers & Pasquini 1997, for A
c type supergiants. Their predicted rate is too low for G and possibly F
c F type giants, however.
c ------------------------------------------------------------------------------
c Input:
c     T:    Effective temperature, in Kelvin
c     logL: log10 L, in solar units
c     logZ: Metallicity
c Returns:
c     mass loss rate, in solar masses per year
c Notes:
c     based on their expression (15). Might be better to use (18) instead,
c     which requires the force parameters alpha and k (Abbott 1982).
c     The metallicity dependence was taken from Vink et al (2001), for O
c     and B star winds. The physical mechanism (line driving) is the same.
c ------------------------------------------------------------------------------
      FUNCTION CALC_MDOT_ACHMAD(T, logL, logZ)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: T, logL, logZ
      DOUBLE PRECISION :: CALC_MDOT_ACHMAD
      DOUBLE PRECISION :: F(10) = (/ -17.359, -16.925, -16.648, -16.255,
     &                               -15.755, -15.241, -14.883, -14.620, 
     &                               -14.430, -14.430 /)
      INTEGER :: I
      DOUBLE PRECISION :: X, LOG_MDOT

      call clear_errno()

! NB: Maybe this mass loss rate should not be trusted below 8000K...
!      IF (T < 5500.0) THEN
!         I = 1
!         X = 0.0
!         call set_errno(err_mdot_extrapolated_low_T)
!         CALC_MDOT_ACHMAD = 0.0d0
!         return
      IF (T < 8000.0) THEN
         I = 5
         X = 0.0
         call set_errno(err_mdot_extrapolated_low_T)
      ELSEIF (T > 9500.0) THEN
         I = 9
         X = 0.0
         call set_errno(err_mdot_extrapolated_high_T)
      ELSE
         I = INT((T-5500.0D0)/500.0D0) + 1
         X = MOD(T, 500.0D0)/500.0
      END IF

      LOG_MDOT = X*F(I+1) + (1.0-X)*F(I)  + 1.646*logL + 0.85*logZ
      CALC_MDOT_ACHMAD = 10**LOG_MDOT
      END FUNCTION

c ------------------------------------------------------------------------------
c CALC_MDOT_VINK:
c Vink mass loss recipe for massive stars including metallicity dependence.
c Based on Vink et al (99, 00, 01) see http://www.astro.keele.ac.uk/~jsv/
c Implemented in STARS: SdM July 2006
c This version: EG September 2008
c ------------------------------------------------------------------------------
c Input parameters:
c     Teff: effective temperature
c     lum:  luminosity
c     Mass: stellar mass
c     logZ: metallicity (0.02 for solar) taken from settings
c Returns:
c     mass loss rate, in solar masses per year
c Notes:
c     for temperatures below the second bistability jump, sets the
c     mdot_error flag to indicate unreliable extrapolation to low Teff
c ------------------------------------------------------------------------------
      FUNCTION CALC_MDOT_VINK(Mass, Lum, Teff, logZ)
      USE CONSTANTS
      IMPLICIT NONE  
      DOUBLE PRECISION :: CALC_MDOT_VINK
! constants   
      DOUBLE PRECISION, PARAMETER :: SMALL = 1.0D-40   
      DOUBLE PRECISION, PARAMETER :: SIGMAE = 0.325
! input/output variabels
      DOUBLE PRECISION, INTENT(IN) :: LUM, MASS, TEFF, logZ
! Coefficients below and above the first bistability jump
      DOUBLE PRECISION :: COEFF_BELOW(7) = (/-6.688, 2.210, -1.339, -1.601, 1.07, 0., 0.85/)
      DOUBLE PRECISION :: COEFF_ABOVE(7) = (/-6.697, 2.194, -1.313, -1.226, 0.933, -10.92, 0.85/)
      DOUBLE PRECISION, PARAMETER :: TEFF_SMOOTH = 1000.0
! parameters
      DOUBLE PRECISION :: COEFF(7), RATIO !=vesc/vinf 
! derived variables
      DOUBLE PRECISION :: GAMMA
      DOUBLE PRECISION :: CHAR_RHO, TEFF_JUMP1, TEFF_JUMP2
      DOUBLE PRECISION :: LOG_TEFF, lOG_MDOT, X
! Use second bistability jump?
      LOGICAL, PARAMETER :: USE_SECOND_JUMP = .TRUE.
! Exponential fall-off rate for cool side of bi-stability jump; motivated
! the Wien approximation for the radiation intensity
      DOUBLE PRECISION :: MDOT_MODULATE

      call clear_errno()

! TODO: SIGMAE should be calculated as in Lamers&Leitherer 1993 instead of being
! constant, although this is probably a minor correction.
      GAMMA = 7.66D-5 *SIGMAE * LUM/MASS
      MDOT_MODULATE = 1.0d0
     
c Determine postions of bistability jumps
      CHAR_RHO   =  -14.94 + ( 3.1875 * Gamma   ) + (0.85 * logZ) !Eq 23, Vink (2001)
      TEFF_JUMP1 = ( 61.2  +   2.59  * CHAR_RHO ) *  1.0D3         !Eq 15, Vink (2001)
      TEFF_JUMP2 = ( 1.0D2 +   6.0D0 * CHAR_RHO ) *  1.0D3         !Eq  6, Vink (2001)

      IF (GAMMA > 0.5) call set_errno(err_mdot_extrapolated_high_L)
      IF (TEFF > 50000.0) call set_errno(err_mdot_extrapolated_high_T)
c Determine postion with respect to temperature jumps 
      IF(TEFF < TEFF_JUMP1) THEN
         RATIO = 1.3
         COEFF(:) = COEFF_BELOW(:)
         IF (TEFF < TEFF_JUMP2) THEN ! below both jumps
! We are probably out of the validity range of this recipe now
            call set_errno(err_mdot_extrapolated_low_T)
            MDOT_MODULATE = DEXP(-TEFF_JUMP2 / TEFF)

! Smooth out the size of the (second) bistability jump over ~ 1000 K
            X = 0.0D0
            IF ( abs(TEFF - TEFF_JUMP2)<TEFF_SMOOTH ) THEN
               X = ( TEFF_JUMP2+TEFF_SMOOTH - TEFF ) / (2.0*TEFF_SMOOTH)
               !X = 0.5*(1.0D0 - COS(X*CPI))   ! Cosine interpolation, smoother
            END IF
            IF (TEFF < TEFF_JUMP2-TEFF_SMOOTH) X = 1
            IF (USE_SECOND_JUMP) COEFF(1) = (1.0d0-X)*COEFF_BELOW(1)+X*(-5.990)
         ENDIF
      ELSE                          ! above both jumps
         ratio = 2.6
         COEFF(:) = COEFF_ABOVE(:)
      ENDIF

! Smooth out the size of the (first) bistability jump over ~ 1000 K
      IF ( abs(TEFF - TEFF_JUMP1)<TEFF_SMOOTH ) THEN
         X = ( TEFF_JUMP1+TEFF_SMOOTH - TEFF ) / (2.0*TEFF_SMOOTH)
         ! X = 0.5*(1.0D0 - COS(X*CPI))   ! Cosine interpolation, smoother
         RATIO = X*1.3 + (1.0D0 - X)*2.6
         COEFF(:) = X*COEFF_BELOW(:) + (1.0D0 - X)*COEFF_ABOVE(:)
      END IF

c get mass loss
      LOG_TEFF = LOG10( MAX(SMALL, TEFF  / 4.0D4) )
      LOG_MDOT = COEFF(1)
     &        + COEFF(2) *  LOG10( MAX(SMALL, LUM   / 1.0D5) )
     &        + COEFF(3) *  LOG10( MAX(SMALL, MASS  / 3.0D1) )
     &        + COEFF(4) *  LOG10( MAX(SMALL, RATIO / 2.0D0) )
     &        + COEFF(5) *  LOG_TEFF
     &        + COEFF(6) *  LOG_TEFF**2
     &        + COEFF(7) *  logZ

      IF (LOG10(LUM) < 4.0) THEN
         LOG_MDOT = LOG_MDOT - 2.0*(4.0 - LOG10(LUM))
         call set_errno(err_mdot_extrapolated_low_L)
      END IF 

      IF (TEFF < 8000.0) THEN
         LOG_MDOT = LOG_MDOT - (8000.0-TEFF)/200.0
      END IF
    
      CALC_MDOT_VINK = 10**lOG_MDOT * MDOT_MODULATE
      END FUNCTION



c ------------------------------------------------------------------------------
! CALC_MDOT_KUDRITZKI:
!  mass loss rate according to Kudritzki 2002 ApJ, based on equations
! (58)-(60) in the paper.
c ------------------------------------------------------------------------------
! Input:
!  logL - log10 luminosity (Lsun)
!  TEFF - Effective temperature (K)
! Returns:
!  Mass loss rate, in solar masses per year
c ------------------------------------------------------------------------------
      FUNCTION CALC_MDOT_KUDRITZKI(logL, TEFF, logZ)
      USE CONSTANTS
      IMPLICIT NONE  
      DOUBLE PRECISION :: CALC_MDOT_KUDRITZKI
      DOUBLE PRECISION, INTENT(IN) :: logL, TEFF, logZ
      ! Temperature bounds for the table
      DOUBLE PRECISION :: TEFF_BOUNDS(3) = (/ 60E3, 50E3, 40E3 /)
      ! Interpolation coeeficients and table data
      DOUBLE PRECISION :: ZMIN_TABLE(3,3) = RESHAPE (
     &         (/  -3.40D0, -0.40D0, -0.65D0, 
     &             -3.85D0, -0.05D0, -0.60D0, 
     &             -4.45D0,  0.35D0, -0.80D0 /), (/3,3/) )
      DOUBLE PRECISION :: QMIN_TABLE(3,3) = RESHAPE (
     &         (/  -8.00D0, -1.20D0, 2.15D0, 
     &            -10.35D0,  3.25D0, 0.00D0, 
     &            -11.75D0,  3.65D0, 0.00D0 /), (/3,3/) )
      DOUBLE PRECISION :: QNIL_TABLE(3,3) = RESHAPE(
     &         (/  -5.99D0,  1.00D0, 1.50D0, 
     &             -4.84D0,  0.50D0, 1.00D0, 
     &             -5.20D0,  0.93D0, 0.85D0 /), (/3,3/) )
      DOUBLE PRECISION :: AZ(3), AQMIN(3), AQNIL(3)
      DOUBLE PRECISION :: TI, TIM
      DOUBLE PRECISION :: L, ZMIN, QMIN, QNIL, LOG_MDOT
      INTEGER :: I
      
      call clear_errno()

      IF (TEFF < TEFF_BOUNDS(3)) MDOT_ERRNO = ERR_MDOT_EXTRAPOLATED_LOW_T
      IF (TEFF > TEFF_BOUNDS(1)) MDOT_ERRNO = ERR_MDOT_EXTRAPOLATED_HIGH_T
      ! Determine coefficients by interpolation, clip results to the edge
      ! of the table.
      IF (TEFF > TEFF_BOUNDS(1)) THEN
         AZ(:) = ZMIN_TABLE(:, 1)
         AQMIN(:) = QMIN_TABLE(:, 1)
         AQNIL(:) = QNIL_TABLE(:, 1)
      ELSEIF (TEFF < TEFF_BOUNDS(3)) THEN
         AZ(:) = ZMIN_TABLE(:, 3)
         AQMIN(:) = QMIN_TABLE(:, 3)
         AQNIL(:) = QNIL_TABLE(:, 3)
      ELSE
         ! Interpolate, determine interpolation factor
         IF (TEFF > TEFF_BOUNDS(2)) THEN
            I = 2
         ELSE
            I = 3
         ENDIF
         TI = (TEFF - TEFF_BOUNDS(I)) / (TEFF_BOUNDS(I-1) - TEFF_BOUNDS(I))
         TIM = 1.0-TI
         ! Interpolate tables between these temperatures
         AZ(:) = TI * ZMIN_TABLE(:, I-1) + TIM * ZMIN_TABLE(:, I)
         AQMIN(:) = TI * QMIN_TABLE(:, I-1) + TIM * QMIN_TABLE(:, I)
         AQNIL(:) = TI * QNIL_TABLE(:, I-1) + TIM * QNIL_TABLE(:, I)
      END IF
      ! Determine ZMIN, QMIN and QNIL as functions of L
      L = logL - 6.0
      ZMIN = AZ(1) + AZ(2)*L + AZ(3) * L**2
      QMIN = AQMIN(1) + AQMIN(2)*L + AQMIN(3) * L**2
      QNIL = AQNIL(1) + AQNIL(2)*L + AQNIL(3) * L**2
      
      LOG_MDOT = (QNIL-QMIN) / SQRT(-ZMIN) * SQRT(logZ - ZMIN) + QMIN
      CALC_MDOT_KUDRITZKI = 10**lOG_MDOT
      END FUNCTION



c ------------------------------------------------------------------------------
! CALC_MDOT_WR_NL:
!  Wolf-Rayet star mass loss rate from Nugis & Lamers 2000 A&A,
!  based on their expressions (20) and (21). An alternative would be their (22)
!  for all WR types, but Nugis & Lamers give a different scaling with initial
!  Z for the different types.
!  The distinction between different Wolf-Rayet types is based on
!  Eldirdge & Vink 2006 A&A
c ------------------------------------------------------------------------------
! Input:
!  LOGL  - log10 luminosity, in solar units
!  LTEFF - log10 effective temperature, in Kelvin
!  X     - Surface abundance of hydrogen, by mass
!  Y     - Surface abundance of helium, by mass
!  XO    - Surface abundance fraction of oxygen, by mass
!  XC    - Surface abundance fraction of carbon, by mass
! Returns:
!  CALC_MDOT_WR_NL - mass loss rate, in solar masses per year
c ------------------------------------------------------------------------------
      FUNCTION CALC_MDOT_WR_NL(LOGL, LTEFF, X, Y, XO, XC, logZ)
      IMPLICIT NONE
      DOUBLE PRECISION :: CALC_MDOT_WR_NL
      DOUBLE PRECISION, INTENT(IN) :: LOGL, LTEFF, X, Y, XO, XC, logZ
      DOUBLE PRECISION, PARAMETER :: X_WNE = 0.05  ! From Nugis & Lamers
      DOUBLE PRECISION, PARAMETER :: Z_WC = 0.03   ! From Eldridge & Vink
      DOUBLE PRECISION, PARAMETER :: Z_WO = 1.0    ! From Eldridge & Vink
      DOUBLE PRECISION :: LOG_MDOT
      DOUBLE PRECISION :: ZETA
!      DOUBLE PRECISION :: CH2(4), CHI(26,9), COM(27), CAN(9), CBN(9)
!      INTEGER :: KZN(9)
!      COMMON /ATDATA/ CH2, CHI, COM, CAN, CBN, KZN
      
      call clear_errno()

      CALC_MDOT_WR_NL = 0.0D0
!     Check that the star is really a Wolf-Rayer star, following Eldridge & Vink
      IF (LTEFF < 4.0D0 .OR. X>0.4D0) then
         mdot_errno = err_mdot_extrapolated_low_T
         RETURN
      END IF

!     Determine Wolf-Rayet subtype. Note that the mass loss rate is the same
!     for WNL and WNE stars, and for WC and WO stars
!     Convert mass fractions to number fractions
!      ZETA = CAN(2)*(XC/CAN(3) + XO/CAN(5)) / Y
      ZETA = 4.0026*(XC/12.0 + XO/14.003) / Y
      IF (X > X_WNE) THEN           ! WNL
         LOG_MDOT = -13.60 + 1.63*LOGL + 2.22*LOG10(Y)
      ELSEIF (ZETA < Z_WC) THEN     ! WNE
         LOG_MDOT = -13.60 + 1.63*LOGL + 2.22*LOG10(Y)
      ELSEIF (ZETA < Z_WO) THEN     ! WC
         LOG_MDOT = -8.30 + 0.84*LOGL + 2.04*LOG10(Y) + 1.04*LOG10(1.0-X-Y)
      ELSE                          ! WO
         LOG_MDOT = -8.30 + 0.84*LOGL + 2.04*LOG10(Y) + 1.04*LOG10(1.0-X-Y)
      END IF
!     Metallicity scaling for WN stars taken from Vink&de Koter 2005
      LOG_MDOT = LOG_MDOT + 0.86 * MAX(logZ, -3.0)
!     If not quite WR, modulate the mass loss rate (so it is turned on smoothly
!     as a star evolves to become a WR star)
      IF (LTEFF < 4.0D0 .OR. X>0.4D0) THEN   ! Not a WR star
         LOG_MDOT = LOG_MDOT - (MAX(X, 0.4)-0.4)*20.0 - 6.0*(4.00-MIN(4.0,LTEFF))
         MDOT_ERRNO = ERR_MDOT_EXTRAPOLATED_LOW_T
      END IF

      CALC_MDOT_WR_NL = 10**lOG_MDOT
      END FUNCTION

      end module
