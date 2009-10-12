      module extra_elements
      
      ! Nitrogen 14 parameters; these are not actually used anymore though...
      ! Well, EN14 is.
      integer, parameter :: NN14 = 16     ! Number of the independent variable
      integer, parameter :: EN14 = 13     ! Number of the nitrogen equation
      logical :: use_n14_eqn = .true.     ! Can disable N14 eqn
      
      ! Magnesium 24 parameters  these are not actually used anymore though...
      integer, parameter :: NMG24 = 21    ! Number of the independent variable
      integer, parameter :: EMG24 = 14    ! Number of magnesium equation
      logical :: use_mg24_eqn = .false.   ! Can disable Mg24 eqn

      ! Equation number for all chemical elements summed
      integer, parameter :: ESUMX = 15    ! Sum of chemical components = 1

      ! Angular momentum transport
      integer, parameter :: EAMT = 42     ! Diffusion equation
      integer, parameter :: FX_RICH = 103 ! Richardson number
      integer, parameter :: FX_DDSI = 104 ! Mixing coefficient for dynamical shear
      integer, parameter :: FX_VES  = 105 ! Eddington-Sweet circulation velocity
      integer, parameter :: FX_VMU  = 106 ! Current from mu gradient (except for grad mu)
      integer, parameter :: FX_HP   = 107 ! Pressure scale-height
      integer, parameter :: FX_DSSI = 108 ! Mixing coefficient for secular shear
      integer, parameter :: FX_SSSI = 109 ! Stability criterion for SSI
      integer, parameter :: FX_SGF  = 114 ! Conversion factor for diffusion coefficients

      ! Total angular momentum, similar to the moment of inertia
      integer, parameter :: ETAM = 43     ! Total angular momentum
      integer, parameter :: NTAM = 42     ! Total angular momentum
      integer, parameter :: FX_AM  = 111  ! Total angular momentum
      integer, parameter :: FX_AMK = 112  ! Total angular momentum gradient
      
      ! domega/dk, differential rotation
      integer, parameter :: FX_RIN  = 113 ! *real* Richardson number
      integer, parameter :: FX_DES  = 122 ! Mixing coefficient for E-S circulation
      integer, parameter :: FX_DGMU = 123 ! Prefactor for grad mu

      ! Gravitational settling
      integer, parameter :: FX_FH  = 115  ! Advective flux hydrogen
      integer, parameter :: FX_FHe = 116  ! Advective flux helium
      integer, parameter :: FX_FC  = 117  ! Advective flux carbon
      integer, parameter :: FX_FN  = 118  ! Advective flux nitrogen
      integer, parameter :: FX_FO  = 119  ! Advective flux ocygen
      integer, parameter :: FX_FNe = 120  ! Advective flux neon
      integer, parameter :: FX_FMg = 121  ! Advective flux magnesium

      ! Net accretion mass flux
      integer, parameter :: FX_MACC = 110 ! Net mass accretion rate, 1e33g/s

      end module
