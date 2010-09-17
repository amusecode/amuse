module massloss
   use real_kind
   use constants
   
   implicit none
   ! Return codes used by the different functions
   ! No problems encountered
   integer, parameter :: err_mdot_ok = 0
   ! Parameters were out of the validity range for this recipe. Proceed with caution!
   integer, parameter :: err_mdot_extrapolated_high_t =  1
   integer, parameter :: err_mdot_extrapolated_low_t  =  2
   integer, parameter :: err_mdot_extrapolated_high_l =  4
   integer, parameter :: err_mdot_extrapolated_low_l  =  8
   integer, parameter :: err_mdot_extrapolated        = 15
   integer, parameter :: err_mdot_not_used            = 16
   
   integer :: mdot_errno
   
   real(double) :: metallicity_scaling_exponent = 0.8 
   
   ! Optionally correct de Jager for clumping (in an approximate way), but
   ! only for higher temperatures.
   real(double) :: de_jager_clump_factor = 1.0
   
   ! Switches for the different mass loss rates
   
   ! General, to be used if none of the others are
   real(double) :: multiplier_dejager = 1.0
   
   ! RGB rates
   real(double) :: multiplier_schroeder = 1.0
   real(double) :: multiplier_reimers = 0.0
   
   ! AGB rates
   real(double) :: multiplier_vasiliadis_wood = 0.0
   real(double) :: multiplier_wachter = 1.0
   
   ! A type super giants
   real(double) :: multiplier_achmad = 0.0
   
   ! Hot stars
   real(double) :: multiplier_vink = 1.0
   
   ! Stars close to the Eddington limit (and WR stars)
   real(double) :: multiplier_kudritzki = 0.0
   real(double) :: multiplier_nl = 1.0
   
   
contains
   
   
   !----------------------------------------------------!
   ! Helper functions for accessing bit fields in errno !
   !----------------------------------------------------!
   
   subroutine clear_errno
      use real_kind
      
      implicit none
      mdot_errno = 0
   end subroutine clear_errno
   
   
   subroutine set_errno(error_flag)
      use real_kind
      
      implicit none
      integer :: error_flag
      mdot_errno = ior(mdot_errno, error_flag)
   end subroutine set_errno
   
   
   function test_flag(flags, bitmask)
      use real_kind
      
      implicit none
      logical :: test_flag
      integer :: flags, bitmask
      test_flag = iand(flags, bitmask) == bitmask
   end function test_flag
   
   
   function test_errno(error_flag)
      use real_kind
      
      implicit none
      integer :: error_flag
      logical :: test_errno
      test_errno = iand(mdot_errno, error_flag) == error_flag
   end function test_errno
   
   
   ! calculate_mdot:
   ! Compute mass loss rate from stellar parameters, using various recipes
   ! and deciding which of those is most appropriate. The transitions
   ! between the various recipes are supposed to be joined together
   ! smoothly.
   ! INPUT:
   !  Teff  - effective temperature, in Kelvin
   !  L     - luminosity, in solar units
   ! (Ledd  - Eddington luminosity, in solar units)
   !  M     - Mass, in solar units
   !  R     - Radius, in solar units
   !  X     - Surface abundance of hydrogen, by mass
   !  Y     - Surface abundance of helium, by mass
   !  XO    - Surface abundance fraction of oxygen, by mass
   !  XC    - Surface abundance fraction of carbon, by mass
   !  Z0    - Initial metallicity (proxy for the iron content)
   ! RETURNS:
   !  The mass loss rate, in solar masses per year
   
   !function calculate_mdot( Teff, l, ledd, m, r, x, y, xo, xc, z0 )
   function calculate_mdot( Teff, l, m, r, x, y, xo, xc, z0 )
      use real_kind
      use constants
      
      implicit none
      !real(double), intent(in) :: Teff, l, ledd, m, r, x, y, xo, xc, z0
      real(double), intent(in) :: Teff, l, m, r, x, y, xo, xc, z0
      real(double) :: calculate_mdot
      
      ! Store results of different mass loss recipes in an array, identify
      ! them using symbolic constants for the array index
      integer, parameter :: dejager = 1
      integer, parameter :: reimers = 2
      integer, parameter :: schroeder = 3
      integer, parameter :: wachter = 4
      integer, parameter :: vasiliadis_wood = 5
      integer, parameter :: vink = 6
      integer, parameter :: kudritzki = 7
      integer, parameter :: achmad = 8
      integer, parameter :: nl = 9
      integer, parameter :: num_recipes = 9
      
      
      ! Local variables and constants
      real(double), parameter :: cmsny = cmsn/csy ! Solar masses per year
      real(double) :: logZ
      real(double) :: logL, logT
      
      ! Outcome of the different mass loss recipes
      real(double) :: mdot_r(num_recipes)
      
      ! Store error return codes
      integer :: recipe_errno(num_recipes)
      
      ! Different masas loss regimes
      real(double) :: mdot_agb, mdot_rgb, mdot_line, mdot_luminous
      integer :: errno_agb, errno_rgb, errno_line, errno_luminous
      real(double) :: mdot_cool, mdot_hot
      integer :: errno_hot, errno_cool
      real(double) :: mdot
      real(double) :: xx
      
      calculate_mdot = 0.0d0
      mdot_r(:) = 0.0d0
      recipe_errno(:) = 0
      logZ = log10(z0/czsn)
      logL = log10(l)
      logT = log10(Teff)
      
      ! Compute contributions to the mass loss from different driving
      ! mechanisms
      
      ! Default prescription if no other applies:
      !  * de Jager et al.
      mdot_errno = err_mdot_not_used
      if (multiplier_dejager>0.0d0)&
           mdot_r(dejager) = multiplier_dejager * calc_mdot_dejager(logL, logT)
      if (test_errno(err_mdot_extrapolated_low_l)) mdot_r(dejager) = 0.0
      recipe_errno(dejager) = mdot_errno
      
      ! Scale the de Jager rate with the *initial* metallicity (ie, Fe
      ! abundance). The other mass loss recipes either have a proper
      ! metallicity dependence (Vink), don't depend on metallicity (Schroeder
      ! & Cuntz) or depend on it in some different way (presumably dust-driven
      ! winds).
      ! Because this mass loss scaling should only apply to line driven winds,
      ! apply it only for Teff>8000, which is roughly where Achmad&al. find
      ! the border of line driven winds.
      if (Teff>9500.0) &
           mdot_r(dejager)=mdot_r(dejager)*(z0/czsn)**metallicity_scaling_exponent
      
      ! Wind clumping correction for the radiation driven portion of the de
      ! Jager rate. Turn on smoothly towards higher T
      if (Teff>8000) then
         xx = min(1.d0, (Teff-8000.d0)/1000.d0)
         mdot_r(dejager) = mdot_r(dejager)*(xx*de_jager_clump_factor + (1.0-xx)*1.0)
      end if
      
      ! Mass loss mechanisms on the giant branch (driven by Alfven waves(?))
      !  * Classical Reimers
      !  * Eggleton version of Reimers
      !  * Schroeder & Cuntz (default)
      mdot_errno = err_mdot_not_used
      if (multiplier_reimers>0.0d0)&
           mdot_r(reimers) = multiplier_reimers * calc_mdot_reimers(Teff, l, m, r)
      recipe_errno(reimers) = mdot_errno
      
      mdot_errno = err_mdot_not_used
      if (multiplier_schroeder>0.0d0)&
           mdot_r(schroeder) = multiplier_schroeder * calc_mdot_schroeder(Teff, l, m, r)
      recipe_errno(schroeder) = mdot_errno
      
      !> \todo FIXME: we don't relly want to include fudges like this here. The
      !! proper thing to do would be to
      !! a. Compute all mass loss rates and store their return flags. Get all
      !! of them to return an "interpolation factor" that is to be used if this
      !! rate needs to be smoothly joined to another recipe. Return flags need
      !! to be a bit field because we need to know if the failure is just in
      !! temperature or luminosity, or both.
      !! b. In each regime, find the best rate. This is the rate that for which
      !! the parameters are "within range". If there are more valid rates, then
      !! take the maximum value of all valid rates. If all of them are
      !! extrapolated, we record that information and again take the maxiumum.
      !! c. If there are no valid rates for this set of parameters we
      !! interpolate.
      !! c2. Determine rates on the "hot" and "cool" side of the current set of
      !! parameters (M, L).
      !! c3. On the cool side: if all available rates fail at high and low
      !! luminosities as well, then interpolate the RGB rate with the AGB rate.
      !! c4. On the hot side: the Vink rate should not fail at low or high
      !! luminosities, although it is forced to 0 at smaller luminosities.
      !! c5. Interpolate between the high T and low T solutions, like this:
      !!  mdot = xh*(1-xl)*mdot_h + (1-xh)*xl*mdotl
      !! c6. If there is a fall-back rate (de Jager):
      !!  return max(mdot_fallback, mdot)
      !<
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
      if (multiplier_wachter>0.0d0)&
           mdot_r(wachter) = multiplier_wachter * calc_mdot_wachter(logT, logL, m)
      recipe_errno(wachter) = mdot_errno
      
      mdot_errno = err_mdot_not_used
      if (multiplier_vasiliadis_wood>0.0d0)&
           mdot_r(vasiliadis_wood) = &
           multiplier_vasiliadis_wood * calc_mdot_vasiliadis_wood(l, m, r)
      recipe_errno(vasiliadis_wood) = mdot_errno
      
      ! Mass loss from A-type supergiants (line driven)
      !  * de Jager et al
      !  * Achmad (default)
      mdot_errno = err_mdot_not_used
      if (multiplier_achmad>0.0d0)&
           mdot_r(achmad) = multiplier_achmad * calc_mdot_achmad(Teff, logL, logZ)
      recipe_errno(achmad) = mdot_errno
      
      ! Mass loss from hot stars (line driven)
      !  * de Jager et al
      !  * Vink (default)
      mdot_errno = err_mdot_not_used
      if (multiplier_vink>0.0d0)&
           mdot_r(vink) = multiplier_vink * calc_mdot_vink(m, l, Teff, logZ)
      recipe_errno(vink) = mdot_errno
      
      ! Mass loss close to the Eddington limit (continuum driven, WR wind)
      !  * Nugis & Lamers 2002 (default)
      !  * de Jager et al.
      !  * Kudritzki 2002
      mdot_errno = err_mdot_not_used
      if (multiplier_kudritzki>0.0d0)&
           mdot_r(kudritzki) = multiplier_kudritzki * calc_mdot_kudritzki(logL, Teff, logZ)
      recipe_errno(kudritzki) = mdot_errno
      
      mdot_errno = err_mdot_not_used
      if (multiplier_nl>0.0d0)&
           mdot_r(nl) = multiplier_nl * calc_mdot_wr_nl(logL, logT, x, y, xo, xc, logZ)
      recipe_errno(nl) = mdot_errno
      
      
      
      ! Determine mass loss rate by picking from the different prescriptions above:
      !> \todo FIXME: this currently fails if both rates are valid (have no errors)
      mdot_agb = max(mdot_r(wachter), mdot_r(vasiliadis_wood))
      errno_agb = recipe_errno(wachter)
      
      mdot_rgb = max(mdot_r(reimers), mdot_r(schroeder))
      errno_rgb = recipe_errno(schroeder)
      if (errno_rgb == err_mdot_not_used) errno_rgb = mdot_r(reimers)
      
      mdot_line = max(mdot_r(vink), mdot_r(achmad))
      errno_line = recipe_errno(vink)
      
      mdot_luminous = 0.0
      errno_luminous = err_mdot_not_used
      !> \todo FIXME: this commented-out bit of code blatently doesn't do what it's
      !! supposed to do... which is decide on which of the two rates for
      !! luminous stars it should use.
      !<
      !      if (recipe_errno(KUDRITZKI) == err_mdot_ok) then
      !         errno_luminous = err_mdot_ok
      !         mdot = max(mdot_r(KUDRITZKI), mdot_luminous)
      !         if (mdot_r(NL) > mdot_luminous .and. recipe_errno(NL) == err_mdot_ok) then
      !            mdot = max(mdot_r(NL), mdot_luminous)
      !         end if
      !      else if (recipe_errno(NL) == err_mdot_ok) then
      !         mdot_luminous = mdot_r(NL)
      !         mdot_luminous = max(mdot_r(NL), mdot_luminous)
      !      else
      !         mdot_luminous = max(mdot_r(KUDRITZKI), mdot_r(NL))
      !         errno_luminous = ior(recipe_errno(KUDRITZKI), recipe_errno(NL))
      !         errno_luminous = iand(errno_luminous, not(err_mdot_not_used))
      !      end if
      mdot_luminous = mdot_r(nl)
      errno_luminous = recipe_errno(nl)
      
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
      if (test_flag(errno_agb, err_mdot_extrapolated_low_l) .and.&
           test_flag(errno_rgb, err_mdot_extrapolated_high_l)) then
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
      if (errno_line == err_mdot_ok .or. &
           test_flag(errno_line, err_mdot_extrapolated_low_l)) then
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
      if (test_flag(errno_hot, err_mdot_extrapolated_low_t) .and.&
           test_flag(errno_cool, err_mdot_extrapolated_high_t) ) then
         !> \todo FIXME: Add interpolation information
         mdot = 0.5*(mdot_cool+mdot_hot)
         mdot_errno = err_mdot_extrapolated
      end if
      if (mdot_errno == err_mdot_not_used) then
         mdot = mdot_cool
         mdot_errno = errno_hot
      end if
      
      if (mdot_errno /= err_mdot_ok) then
         mdot = max(mdot_r(dejager), mdot)
      end if
      !mdot = max(mdot_agb, mdot_rgb, mdot_line, mdot_luminous)
      
      ! If the adopted rate turns out to be NL but that was extrapolated, max
      ! it with the de Jager rate.
      !      if (test_errno(err_mdot_extrapolated_low_T) .and. mdot==mdot_r(NL))
      !     &   mdot = max(mdot, mdot_r(DEJAGER))
      
      ! Fall back to the default rate in case none of those above apply
      !      if (mdot < 1.0d-10) mdot = max(mdot, mdot_r(DEJAGER))
      
      !> \todo FIXME: for hot, bright low mass stars, force a high mass loss rate.
      !! This is to make sure the envelope of an AGB star is stripped quickly
      !! once its on its way to become a white dwarf.
      !!if (M<1.0 .and. logL>3.7 .and. Teff>4000.0) mdot = 1.0e-4
      !<
      calculate_mdot = mdot
   end function calculate_mdot
   
   
   
   ! ------------------------------------------------------------------------------
   ! CALC_MDOT_REIMERS:
   !  Calculate the mass loss rate for red giants, using the mass loss
   !  prescription from Reimers (1975)
   ! ------------------------------------------------------------------------------
   ! Input parameters:
   !     T: Effective temperature, in K (unused, except to check validity range)
   !     L: Luminosity, in solar units
   !     M: Mass, in solar units
   !     R: Radius, in solar units
   ! Returns:
   !     Mass loss rate, in solar masses per year
   ! ------------------------------------------------------------------------------
   function calc_mdot_reimers ( t, l, m, r )
      use real_kind
      use constants
      
      implicit none
      real(double), intent(in) :: l, m, r, t
      real(double) :: calc_mdot_reimers
      real(double) :: lrm, tcrit, suppress
      
      call clear_errno()
      
      ! At higher effecive temperatures (above say 5000K) this mass loss rate
      ! no longer applies. Turn it off smoothly by modulating it with a term
      ! that suppresses it when T is large.
      ! By itself the mass loss rate would diverge.
      tcrit = min(5000.d0, 5000.d0 - log10(l)*400.d0)
      suppress = 1.0
      if (t>tcrit) then
         !suppress = -1.0/(Tcrit - T)
         !call set_errno(err_mdot_extrapolated_high_T)
      end if
      
      lrm = l*r/m * suppress
      
      calc_mdot_reimers = 4.0e-13*lrm
   end function calc_mdot_reimers
   
   
   
   ! ------------------------------------------------------------------------------
   ! CALC_CRITICAL_RGB_LUMINOSITY
   ! Determine minimum critical luminosity for driving a dust-driven wind, see
   ! Schroeder & al (1999).
   ! Fortunately, this diverges towards high temperatures, correctly
   ! stating that very hot stars don't drive a dust-driven wind.
   !> \todo FIXME: should use a constraint involving the dust opacity instead, see
   !! Wachter & al.
   !<
   ! ------------------------------------------------------------------------------
   ! Input parameters:
   !     T: Effective temperature, in Kelvin
   !     logM: logarithm of the Mass, in solar units
   ! Returns:
   !     Critical luminosity, in solar units
   ! ------------------------------------------------------------------------------
   function calc_critical_rgb_luminosity ( t, logM, m )
      use real_kind
      use constants
      
      implicit none
      real(double), intent(in) :: t, logM, m
      real(double) :: calc_critical_rgb_luminosity
      real(double) :: lcrit, x
      
      if (t<2800.0) then
         lcrit = log10(((t-2600.0)*5.0 + 5000.0)) + logM
      else
         x = (t-2800.0)/200.0
         lcrit = log10( (1.0-x)*6000.0 + x*1.0d4*(m-0.05)) + logM
      end if
      calc_critical_rgb_luminosity = lcrit
   end function calc_critical_rgb_luminosity
   
   
   
   ! ------------------------------------------------------------------------------
   ! CALC_MDOT_SCHROEDER:
   !  Calculate the mass loss rate for red giants, using the mass loss
   !  prescription from Schroeder&Cuntz (2005).
   !  Mass loss is driven by Alfven waves and does NOT depend on metallicity
   ! ------------------------------------------------------------------------------
   ! Input parameters:
   !     T: Effective temperature, in Kelvin
   !     L: Luminosity, in solar units
   !     M: Mass, in solar units
   !     R: Radius, in solar units
   ! Returns:
   !     Mass loss rate, in solar masses per year
   ! ------------------------------------------------------------------------------
   function calc_mdot_schroeder ( t, l, m, r )
      use real_kind
      use constants
      
      implicit none
      real(double), intent(in) :: t, l, m, r
      real(double) :: calc_mdot_schroeder
      real(double) :: g_ratio, lrm, lcrit, tcrit, suppress
      
      call clear_errno()
      
      ! At higher effecive temperatures (above say 5000K) this mass loss rate
      ! no longer applies. Turn it off smoothly by modulating it with a term
      ! that suppresses it when T is large.
      ! By itself the mass loss rate would diverge.
      tcrit = min(5000.d0, 5000.d0 - log10(l)*400.d0)
      !Tcrit = 5000.0 - log10(L)*400.0
      suppress = 1.0
      if (t>tcrit) then
         suppress = -100.0/(tcrit - t)
         call set_errno(err_mdot_extrapolated_high_t)
      end if
      
      lcrit = calc_critical_rgb_luminosity(t, log10(m), m)
      if (log10(l) > 1.2*lcrit) then
         call set_errno(err_mdot_extrapolated_high_l)
      end if
      
      ! Surface gravity, in solar units (ignores centrifugal term!)
      g_ratio = r**2/m
      lrm = l*r/m * suppress
      
      calc_mdot_schroeder = 8.0d-14*lrm*(t/4.0d3)**3.5 * (1.0d0+g_ratio/(4.3d3))
   end function calc_mdot_schroeder
   
   
   
   ! ------------------------------------------------------------------------------
   ! CALC_MDOT_WACHTER:
   !  Calculate the mass loss rate for AGB stars, from Wachter & al. (2002)
   !  Mass loss is a dust-driven superwind
   ! ------------------------------------------------------------------------------
   ! Input parameters:
   !     logT: log10 effective temperature, in Kelvin
   !     logL: log10 luminosity, in solar units
   !     M: Mass, in solar units
   ! Returns:
   !     Mass loss rate, in solar masses per year
   ! ------------------------------------------------------------------------------
   function calc_mdot_wachter ( logT, logL, m )
      use real_kind
      use constants
      
      implicit none
      real(double), intent(in) :: logT, logL, m
      real(double) :: calc_mdot_wachter
      real(double) :: log_mdot, lcrit, t, logM, suppress
      
      call clear_errno()
      
      ! Determine minimum critical luminosity for driving a dust-driven wind, see
      ! Schroeder & al (1999).
      ! Fortunately, this diverges towards high temperatures, correctly
      ! stating that very hot stars don't drive a dust-driven wind.
      suppress = 0.0d0
      t = 10**logT
      logM = log10(m)
      lcrit = calc_critical_rgb_luminosity(t, logM, m)
      if (logL < lcrit) then
         call set_errno(err_mdot_extrapolated_low_l)
         suppress = -5.0*(logL-lcrit)
      end if
      
      ! Mass loss rate based on equation (2) in the paper, where the pulsation
      ! period has been eliminated.
      log_mdot = -4.52 -6.81*(logT-3.41) + 2.47*(logL-4) -1.95*logM -suppress
      calc_mdot_wachter = 10**log_mdot
   end function calc_mdot_wachter
   
   
   
   ! ------------------------------------------------------------------------------
   ! CALC_MDOT_VASILIADIS_WOOD
   !  Calculate the mass loss rate for AGB stars, from Vasiliadis&Wood 1993
   !  Mass loss is a pulsation driven superwind
   ! ------------------------------------------------------------------------------
   ! Input parameters:
   !     L: log10 luminosity, in solar units
   !     M: Mass, in solar units
   !     R: Radius, in solar units
   ! Returns:
   !     Mass loss rate, in solar masses per year
   ! ------------------------------------------------------------------------------
   function calc_mdot_vasiliadis_wood ( l, m, r )
      use real_kind
      use constants
      
      implicit none
      real(double), intent(in) :: l, m, r
      real(double) :: calc_mdot_vasiliadis_wood
      real(double) :: pp, max_mdot, log_mdot, vw
      real(double) :: tcrit, t
      
      call clear_errno()
      
      ! Maximum value of the mass loss rate, from equation (1)
      max_mdot = 1.36e-9*l
      
      ! Find Mira pulsation period, equation (4)
      pp = 10**(-2.07 + 1.94*r - 0.9*m)
      
      ! Wind velocity, equation (3), in km/s
      vw = -13.5 + 0.056*pp
      if (vw < 3.0d0 .or. l < 1.0e3) then
         call set_errno(err_mdot_extrapolated_low_l)
      end if
      
      ! Find mass loss rate, from equation (2) (equation 5 for masses above
      ! 2.5)
      log_mdot = -11.4d0 + 0.0125d0*(pp-100.0d0*max(m-2.5d0,0.0d0))
      
      ! Suppress V&W mass loss rate far from the Hayashi track
      tcrit = min(5000.d0, 5000.d0 - log10(l)*400.d0)
      !Tcrit = 5000.0 - log10(L)*400.0
      t = (l/r**2)**0.25 * 5770
      if (t>tcrit) then
         call set_errno(err_mdot_extrapolated_high_t)
      end if
      
      ! mdot is given by the minimum of eq. (1) and eq. (5)
      calc_mdot_vasiliadis_wood = min(max_mdot, 10**log_mdot)
   end function calc_mdot_vasiliadis_wood
   
   
   
   ! ------------------------------------------------------------------------------
   ! CALC_MDOT_DEJAGER:
   !  Empirical mass loss rate for luminous star, de Jager et al (1988)
   !  This mass loss rate should scale with the metallicity
   ! ------------------------------------------------------------------------------
   ! Input parameters:
   !     TT: log10 effective temperature/K
   !     LL: luminosity
   ! Returns:
   !     Mass loss rate, in solar masses per year
   ! ------------------------------------------------------------------------------
   function calc_mdot_dejager ( ll, tt )
      use real_kind
      
      implicit none
      real(double) :: calc_mdot_dejager,ll,tt
      
      integer :: i,j
      real(double) :: mdot
      
      real(double) :: a(6,6)
      data a/ &
           6.34916d0,3.41678d0,-1.08683d0, 0.13095d0,0.22427d0,0.11968d0,  &!J=1
           -5.04240d0,0.15629d0, 0.41952d0,-0.09825d0,0.46591d0,0.0d0,     &!J=2
           -0.83426d0,2.96244d0,-1.37272d0, 0.13025d0,2*0.0d0,             &!J=3
           -1.13925d0,0.33659d0,-1.07493d0, 3*0.0d0,                       &!J=4
           -0.12202d0,0.57576d0, 10*0.0d0/                                  !J=5,6
      
      real(double) :: cheb
      
      
      call clear_errno()
      
      mdot = 0.0d0
      if ( tt < 3.301d0 ) mdot_errno = err_mdot_extrapolated_low_t
      if ( tt > 4.799d0 ) mdot_errno = err_mdot_extrapolated_high_t
      if ( ll < 2.501d0 ) mdot_errno = err_mdot_extrapolated_low_l
      if ( ll > 6.699d0 ) mdot_errno = err_mdot_extrapolated_high_l
      do i = 0, 5
         do j = 0, 5 - i
            mdot = mdot + a(i + 1, j + 1)*cheb((tt - 4.05d0)/0.75d0, i)&
                 *cheb((ll - 4.60d0)/2.10d0, j)
         end do
      end do
      calc_mdot_dejager = 10**(-mdot)
   end function calc_mdot_dejager
   
   
   
   
   
   ! ------------------------------------------------------------------------------
   ! CALC_MDOT_ACHMAD:
   ! Calculate the mass loss rate of Achmad, Lamers & Pasquini 1997, for A
   ! type supergiants. Their predicted rate is too low for G and possibly F
   ! F type giants, however.
   ! ------------------------------------------------------------------------------
   ! Input:
   !     T:    Effective temperature, in Kelvin
   !     logL: log10 L, in solar units
   !     logZ: Metallicity
   ! Returns:
   !     mass loss rate, in solar masses per year
   ! Notes:
   !     based on their expression (15). Might be better to use (18) instead,
   !     which requires the force parameters alpha and k (Abbott 1982).
   !     The metallicity dependence was taken from Vink et al (2001), for O
   !     and B star winds. The physical mechanism (line driving) is the same.
   ! ------------------------------------------------------------------------------
   function calc_mdot_achmad(t, logL, logZ)
      use real_kind
      
      implicit none
      real(double), intent(in) :: t, logL, logZ
      real(double) :: calc_mdot_achmad
      real(double) :: f(10) = (/ -17.359, -16.925, -16.648, -16.255,&
           -15.755, -15.241, -14.883, -14.620, &
           -14.430, -14.430 /)
      integer :: i
      real(double) :: x, log_mdot
      
      call clear_errno()
      
      ! NB: Maybe this mass loss rate should not be trusted below 8000K...
      !      IF (T < 5500.0) THEN
      !         I = 1
      !         X = 0.0
      !         call set_errno(err_mdot_extrapolated_low_T)
      !         CALC_MDOT_ACHMAD = 0.0d0
      !         return
      if (t < 8000.0) then
         i = 5
         x = 0.0
         call set_errno(err_mdot_extrapolated_low_t)
      else if (t > 9500.0) then
         i = 9
         x = 0.0
         call set_errno(err_mdot_extrapolated_high_t)
      else
         i = int((t-5500.0d0)/500.0d0) + 1
         x = mod(t, 500.0d0)/500.0
      end if
      
      log_mdot = x*f(i+1) + (1.0-x)*f(i)  + 1.646*logL + 0.85*logZ
      calc_mdot_achmad = 10**log_mdot
   end function calc_mdot_achmad
   
   
   
   ! ------------------------------------------------------------------------------
   ! CALC_MDOT_VINK:
   ! Vink mass loss recipe for massive stars including metallicity dependence.
   ! Based on Vink et al (99, 00, 01) see http://www.astro.keele.ac.uk/~jsv/
   ! Implemented in STARS: SdM July 2006
   ! This version: EG September 2008
   ! ------------------------------------------------------------------------------
   ! Input parameters:
   !     Teff: effective temperature
   !     lum:  luminosity
   !     Mass: stellar mass
   !     logZ: metallicity (0.02 for solar) taken from settings
   ! Returns:
   !     mass loss rate, in solar masses per year
   ! Notes:
   !     for temperatures below the second bistability jump, sets the
   !     mdot_error flag to indicate unreliable extrapolation to low Teff
   ! ------------------------------------------------------------------------------
   function calc_mdot_vink(mass, lum, Teff, logZ)
      use real_kind
      use constants
      
      implicit none  
      real(double) :: calc_mdot_vink
      ! constants   
      real(double), parameter :: small = 1.0d-40   
      real(double), parameter :: sigmae = 0.325
      ! input/output variabels
      real(double), intent(in) :: lum, mass, Teff, logZ
      ! Coefficients below and above the first bistability jump
      real(double) :: coeff_below(7) = (/-6.688, 2.210, -1.339, -1.601, 1.07, 0., 0.85/)
      real(double) :: coeff_above(7) = (/-6.697, 2.194, -1.313, -1.226, 0.933, -10.92, 0.85/)
      real(double), parameter :: Teff_smooth = 1000.0
      ! parameters
      real(double) :: coeff(7), ratio !=vesc/vinf 
      ! derived variables
      real(double) :: gamma
      real(double) :: char_rho, Teff_jump1, Teff_jump2
      real(double) :: log_Teff, log_mdot, x
      ! Use second bistability jump?
      logical, parameter :: use_second_jump = .true.
      ! Exponential fall-off rate for cool side of bi-stability jump; motivated
      ! the Wien approximation for the radiation intensity
      real(double) :: mdot_modulate
      
      call clear_errno()
      
      ! TODO: SIGMAE should be calculated as in Lamers&Leitherer 1993 instead of being
      ! constant, although this is probably a minor correction.
      gamma = 7.66d-5 *sigmae * lum/mass
      mdot_modulate = 1.0d0
      
      ! Determine postions of bistability jumps
      char_rho   =  -14.94 + ( 3.1875 * gamma   ) + (0.85 * logZ) !Eq 23, Vink (2001)
      Teff_jump1 = ( 61.2  +   2.59  * char_rho ) *  1.0d3         !Eq 15, Vink (2001)
      Teff_jump2 = ( 1.0d2 +   6.0d0 * char_rho ) *  1.0d3         !Eq  6, Vink (2001)
      
      if (gamma > 0.5) call set_errno(err_mdot_extrapolated_high_l)
      if (Teff > 50000.0) call set_errno(err_mdot_extrapolated_high_t)
      ! Determine postion with respect to temperature jumps 
      if(Teff < Teff_jump1) then
         ratio = 1.3
         coeff(:) = coeff_below(:)
         if (Teff < Teff_jump2) then ! below both jumps
            ! We are probably out of the validity range of this recipe now
            call set_errno(err_mdot_extrapolated_low_t)
            mdot_modulate = dexp(-Teff_jump2 / Teff)
            
            ! Smooth out the size of the (second) bistability jump over ~ 1000 K
            x = 0.0d0
            if ( abs(Teff - Teff_jump2)<Teff_smooth ) then
               x = ( Teff_jump2+Teff_smooth - Teff ) / (2.0*Teff_smooth)
               !X = 0.5*(1.0D0 - COS(X*CPI))   ! Cosine interpolation, smoother
            end if
            if (Teff < Teff_jump2-Teff_smooth) x = 1
            if (use_second_jump) coeff(1) = (1.0d0-x)*coeff_below(1)+x*(-5.990)
         end if
      else                          ! above both jumps
         ratio = 2.6
         coeff(:) = coeff_above(:)
      end if
      
      ! Smooth out the size of the (first) bistability jump over ~ 1000 K
      if ( abs(Teff - Teff_jump1)<Teff_smooth ) then
         x = ( Teff_jump1+Teff_smooth - Teff ) / (2.0*Teff_smooth)
         ! X = 0.5*(1.0D0 - COS(X*CPI))   ! Cosine interpolation, smoother
         ratio = x*1.3 + (1.0d0 - x)*2.6
         coeff(:) = x*coeff_below(:) + (1.0d0 - x)*coeff_above(:)
      end if
      
      ! get mass loss
      log_Teff = log10( max(small, Teff  / 4.0d4) )
      log_mdot = coeff(1)&
           + coeff(2) *  log10( max(small, lum   / 1.0d5) )&
           + coeff(3) *  log10( max(small, mass  / 3.0d1) )&
           + coeff(4) *  log10( max(small, ratio / 2.0d0) )&
           + coeff(5) *  log_Teff&
           + coeff(6) *  log_Teff**2&
           + coeff(7) *  logZ
      
      if (log10(lum) < 4.0) then
         log_mdot = log_mdot - 2.0*(4.0 - log10(lum))
         call set_errno(err_mdot_extrapolated_low_l)
      end if
      
      if (Teff < 8000.0) then
         log_mdot = log_mdot - (8000.0-Teff)/200.0
      end if
      
      calc_mdot_vink = 10**log_mdot * mdot_modulate
   end function calc_mdot_vink
   
   
   
   ! ------------------------------------------------------------------------------
   ! CALC_MDOT_KUDRITZKI:
   !  mass loss rate according to Kudritzki 2002 ApJ, based on equations
   ! (58)-(60) in the paper.
   ! ------------------------------------------------------------------------------
   ! Input:
   !  logL - log10 luminosity (Lsun)
   !  TEFF - Effective temperature (K)
   ! Returns:
   !  Mass loss rate, in solar masses per year
   ! ------------------------------------------------------------------------------
   function calc_mdot_kudritzki(logL, Teff, logZ)
      use real_kind
      use constants
      
      implicit none  
      real(double) :: calc_mdot_kudritzki
      real(double), intent(in) :: logL, Teff, logZ
      ! Temperature bounds for the table
      real(double) :: Teff_bounds(3) = (/ 60e3, 50e3, 40e3 /)
      ! Interpolation coeeficients and table data
      real(double) :: zmin_table(3,3) = reshape (          &
           (/  -3.40d0, -0.40d0, -0.65d0,                  &
           -3.85d0, -0.05d0, -0.60d0,                  &
           -4.45d0,  0.35d0, -0.80d0 /), (/3,3/) )
      real(double) :: qmin_table(3,3) = reshape (          &
           (/  -8.00d0, -1.20d0, 2.15d0,                   &
           -10.35d0,  3.25d0, 0.00d0,                   &
           -11.75d0,  3.65d0, 0.00d0 /), (/3,3/) )
      real(double) :: qnil_table(3,3) = reshape(           &
           (/  -5.99d0,  1.00d0, 1.50d0,                   &
           -4.84d0,  0.50d0, 1.00d0,                   &
           -5.20d0,  0.93d0, 0.85d0 /), (/3,3/) )
      real(double) :: az(3), aqmin(3), aqnil(3)
      real(double) :: ti, tim
      real(double) :: l, zmin, qmin, qnil, log_mdot
      integer :: i
      
      call clear_errno()
      
      if (Teff < Teff_bounds(3)) mdot_errno = err_mdot_extrapolated_low_t
      if (Teff > Teff_bounds(1)) mdot_errno = err_mdot_extrapolated_high_t
      ! Determine coefficients by interpolation, clip results to the edge
      ! of the table.
      if (Teff > Teff_bounds(1)) then
         az(:) = zmin_table(:, 1)
         aqmin(:) = qmin_table(:, 1)
         aqnil(:) = qnil_table(:, 1)
      else if (Teff < Teff_bounds(3)) then
         az(:) = zmin_table(:, 3)
         aqmin(:) = qmin_table(:, 3)
         aqnil(:) = qnil_table(:, 3)
      else
         ! Interpolate, determine interpolation factor
         if (Teff > Teff_bounds(2)) then
            i = 2
         else
            i = 3
         end if
         ti = (Teff - Teff_bounds(i)) / (Teff_bounds(i-1) - Teff_bounds(i))
         tim = 1.0-ti
         ! Interpolate tables between these temperatures
         az(:) = ti * zmin_table(:, i-1) + tim * zmin_table(:, i)
         aqmin(:) = ti * qmin_table(:, i-1) + tim * qmin_table(:, i)
         aqnil(:) = ti * qnil_table(:, i-1) + tim * qnil_table(:, i)
      end if
      ! Determine ZMIN, QMIN and QNIL as functions of L
      l = logL - 6.0
      zmin = az(1) + az(2)*l + az(3) * l**2
      qmin = aqmin(1) + aqmin(2)*l + aqmin(3) * l**2
      qnil = aqnil(1) + aqnil(2)*l + aqnil(3) * l**2
      
      log_mdot = (qnil-qmin) / sqrt(-zmin) * sqrt(logZ - zmin) + qmin
      calc_mdot_kudritzki = 10**log_mdot
   end function calc_mdot_kudritzki
   
   
   
   ! ------------------------------------------------------------------------------
   ! CALC_MDOT_WR_NL:
   !  Wolf-Rayet star mass loss rate from Nugis & Lamers 2000 A&A,
   !  based on their expressions (20) and (21). An alternative would be their (22)
   !  for all WR types, but Eldridge & Lamers give a different scaling with initial
   !  Z for the different types.
   !  The distinction between different Wolf-Rayet types is based on
   !  Eldirdge & Vink 2006 A&A
   ! ------------------------------------------------------------------------------
   ! Input:
   !  logL  - log10 luminosity, in solar units
   !  lTeff - log10 effective temperature, in Kelvin
   !  X     - Surface abundance of hydrogen, by mass
   !  Y     - Surface abundance of helium, by mass
   !  XO    - Surface abundance fraction of oxygen, by mass
   !  XC    - Surface abundance fraction of carbon, by mass
   ! Returns:
   !  CALC_MDOT_WR_NL - mass loss rate, in solar masses per year
   ! ------------------------------------------------------------------------------
   function calc_mdot_wr_nl(logL, lTeff, x, y, xo, xc, logZ)
      use real_kind
      
      implicit none
      real(double) :: calc_mdot_wr_nl
      real(double), intent(in) :: logL, lTeff, x, y, xo, xc, logZ
      real(double), parameter :: x_wne = 0.05  ! From Nugis & Lamers
      real(double), parameter :: z_wc = 0.03   ! From Eldridge & Vink
      real(double), parameter :: z_wo = 1.0    ! From Eldridge & Vink
      real(double) :: log_mdot
      real(double) :: zeta
      
      call clear_errno()
      
      calc_mdot_wr_nl = 0.0d0
      !     Check that the star is really a Wolf-Rayer star, following Eldridge & Vink
      if (lTeff < 4.0d0 .or. x>0.4d0) then
         mdot_errno = err_mdot_extrapolated_low_t
         return
      end if
      
      !     Determine Wolf-Rayet subtype. Note that the mass loss rate is the same
      !     for WNL and WNE stars, and for WC and WO stars
      !     Convert mass fractions to number fractions
      !      ZETA = CAN(2)*(XC/CAN(3) + XO/CAN(5)) / Y
      zeta = 4.0026*(xc/12.0 + xo/14.003) / y
      if (x > x_wne) then           ! WNL
         log_mdot = -13.60 + 1.63*logL + 2.22*log10(y)
      else if (zeta < z_wc) then     ! WNE
         log_mdot = -13.60 + 1.63*logL + 2.22*log10(y)
      else if (zeta < z_wo) then     ! WC
         log_mdot = -8.30 + 0.84*logL + 2.04*log10(y) + 1.04*log10(1.0-x-y)
      else                          ! WO
         log_mdot = -8.30 + 0.84*logL + 2.04*log10(y) + 1.04*log10(1.0-x-y)
      end if
      !     Metallicity scaling for WN stars taken from Vink&de Koter 2005
      log_mdot = log_mdot + 0.86d0 * max(logZ, -3.d0)
      !     If not quite WR, modulate the mass loss rate (so it is turned on smoothly
      !     as a star evolves to become a WR star)
      if (lTeff < 4.0d0 .or. x>0.4d0) then   ! Not a WR star
         log_mdot = log_mdot - (max(x, 0.4d0)-0.4d0)*20.d0 - 6.d0*(4.d00-min(4.d0,lTeff))
         mdot_errno = err_mdot_extrapolated_low_t
      end if
      
      calc_mdot_wr_nl = 10**log_mdot
   end function calc_mdot_wr_nl
   
end module massloss



function cheb(x, i)
   use real_kind
   
   implicit none
   real(double) cheb,x
   integer :: i
   
   cheb = cos(i * acos(min(max(x, -1.0_dbl), 1.0_dbl)))
end function cheb


