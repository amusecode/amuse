module extra_elements
   use real_kind
  
   implicit none
   
   ! Flag to detect whether Mg24 is actually being solved for
   logical :: use_mg24_eqn = .false.
  
   ! Equation numbers
   integer, parameter :: EN14 = 13     ! Number of the nitrogen equation
   integer, parameter :: EMG24 = 14    ! Number of magnesium equation
   integer, parameter :: ESUMX = 15    ! Sum of chemical components = 1
   integer, parameter :: EAMT = 42     ! Angular momentum transport
   integer, parameter :: ETAM = 43     ! Total angular momentum
   integer, parameter :: ESi28 = 44    ! Number of silicon equation 
   integer, parameter :: EFe56 = 45    ! Number of iron equation 
  
   ! Index of independent variables in H array
   integer, parameter :: NMG24 = 41    ! Index of Mg24 in H
   integer, parameter :: NSI28 = 42    ! Index of Si28 in H
   integer, parameter :: NFE56 = 43    ! Index of Fe56 in H
   integer, parameter :: NTAM = 44     ! Total angular momentum
  
   ! Extra functions (output from funcs1), ~by category
  
   ! Net accretion mass flux
   integer, parameter :: FX_MACC = 101 ! Net mass accretion rate, 1e33g/s
  
   ! Total angular momentum, similar to the moment of inertia
   integer, parameter :: FX_AM  = 102  ! Total angular momentum
   integer, parameter :: FX_AMK = 103  ! Total angular momentum gradient
   integer, parameter :: FX_SAMK = 104 ! Specific angular momentum
  
   ! domega/dk, differential rotation
   integer, parameter :: FX_RIN  = 105 ! *real* Richardson number
   integer, parameter :: FX_DES  = 106 ! Mixing coefficient for E-S circulation
   integer, parameter :: FX_DGMU = 107 ! Prefactor for grad mu
   ! Angular momentum transport
   integer, parameter :: FX_RICH = 108 ! Richardson number
   integer, parameter :: FX_DDSI = 109 ! Mixing coefficient dynamical shear
   integer, parameter :: FX_VES  = 110 ! Eddington-Sweet circulation velocity
   integer, parameter :: FX_VMU  = 111 ! Current from mu gradient (except for grad mu)
   integer, parameter :: FX_HP   = 112 ! Pressure scale-height
   integer, parameter :: FX_DSSI = 113 ! Mixing coefficient for secular shear
   integer, parameter :: FX_SSSI = 114 ! Stability criterion for SSI
   integer, parameter :: FX_SGF  = 115 ! Conversion factor for diffusion coefficients
  
   ! Gravitational settling
   ! To do this properly, we need 9 fluxes for the 9 main isotopes and
   ! 9x9 = 81 relative diffusion coefficients. For the nucleosynthesis
   ! this is even worse; we'd need a whopping 1640 diffusion terms....
   ! The only way around this is to export the composition gradients to
   ! FUNCS1/FUNCS2, maybe through the semi-implicit functions.
   integer, parameter :: FX_FH  = 116  ! Advective flux hydrogen
   integer, parameter :: FX_FHe = 117  ! Advective flux helium
   integer, parameter :: FX_FC  = 118  ! Advective flux carbon
   integer, parameter :: FX_FN  = 119  ! Advective flux nitrogen
   integer, parameter :: FX_FO  = 120  ! Advective flux ocygen
   integer, parameter :: FX_FNe = 121  ! Advective flux neon
   integer, parameter :: FX_FMg = 122  ! Advective flux magnesium
   integer, parameter :: FX_FSi = 123  ! Advective flux silicon
   integer, parameter :: FX_FFe = 124  ! Advective flux iron

   ! Additional isotopes (Mg, Si, Fe)
   integer, parameter :: FX_X24  = 125  ! Magnesium abundance
   integer, parameter :: FX_X24t = 126  ! Magnesium abundance derivative
   integer, parameter :: FX_X28  = 127  ! Magnesium abundance
   integer, parameter :: FX_X28t = 128  ! Magnesium abundance derivative
   integer, parameter :: FX_X56  = 129  ! Magnesium abundance
   integer, parameter :: FX_X56t = 130  ! Magnesium abundance derivative
  
end module extra_elements
