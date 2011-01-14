! This module contains all mathematical and physical constants
! This does not include constants which are configurable through input files,
module constants
  use real_kind
  
  implicit none
  
  ! Mathematical constants: pi, 4pi, log(10) and 1/3
  real(double), parameter :: CPI = 3.1415926535897932384626433832795029D0
  real(double), parameter :: CPI4 = 4.0D0 * CPI
  real(double), parameter :: CLN = 2.3025850929940456840179914546843642D0
  real(double), parameter :: C3RD = 1.0D0/3.0D0
  
  ! Physical constants
  real(double), parameter :: CL = 2.99792458D10     ! Speed of light
  real(double), parameter :: PLANCK = 6.6260693D-27 ! Planck constant
  real(double), parameter :: CG = 6.6742D-8         ! Newton constant
  real(double), parameter :: BOLTZM = 1.3806505D-16 ! Boltzmann constant
  real(double), parameter :: ECHAR = 4.80320441D-10 ! Electrostatic charge
  
  real(double), parameter :: AMU = 1.660538862D-24  ! Atomic mass unit
  real(double), parameter :: AME = 9.1093826D-28    ! Electron mass
  
  ! Astrophysical parameters
  real(double), parameter :: CMSN = 1.98844D0       ! Solar mass
  real(double), parameter :: CLSN = 3.844D0         ! Solar luminosity
  real(double), parameter :: CRSN = 0.69598D0       ! Solar radius
  real(double), parameter :: CASN = 4.57D9          ! Solar age
  real(double), parameter :: CZSN = 0.02            ! solar metallicity 
  
  ! Units
  real(double), parameter :: CSY = 3.155692597D7    ! Seconds per year
  real(double), parameter :: CSDAY = 86400D0        ! Seconds per day
  real(double), parameter :: EVOLT = ECHAR/CL*1.0D8 ! erg/electron volt
  
  ! Derived constants
  real(double), parameter :: CA = 8.0D0*CPI**5*BOLTZM**4 / (15.0D0 * (CL*PLANCK)**3)
  real(double), parameter :: CME = 1.0D6*EVOLT/AMU  ! (erg/MeV) / (gram/amu)
  real(double), parameter :: CEVB = EVOLT/BOLTZM 
  real(double), parameter :: CR = BOLTZM/AMU        ! gas constant
  
  real(double), parameter :: LAMC = PLANCK/(AME*CL) ! Compton length of electron
  real(double), parameter :: CRHO = 8.0D0*CPI/LAMC**3
  
  real(double), parameter :: CB = CRHO*AME*CL**2 
  real(double), parameter :: CD = CRHO*AMU 
  
  real(double), parameter :: CTE = BOLTZM/(AME*CL**2) 
  real(double), parameter :: CGRT = 6.4D0*(CPI/4.32D4)**6*(1.0D11/CL)**5 
  
  ! Constants that cannot be computed at compile time
  real(double) :: CEN
  real(double) :: CPL
  real(double) :: CG1
  real(double) :: CG2
  
  !Debug:
  integer, parameter :: debug = 0                    ! Debugging info:  0 - no debugging, 2 - print routine names when called
  
contains
  
   
   subroutine initialise_constants
      use real_kind
      use paquette_coefficients
      
      implicit none
      call initialise_collision_integrals
      CEN = (PLANCK/(2.0D0*CPI*AMU*BOLTZM)**0.5D0)**3 / AMU 
      CPL = (CPI4/(AMU*BOLTZM**3))**0.5D0 * ECHAR**3 
      CG1 = 1.0D5*CG**0.5D0 
      CG2 = CG1*(CSDAY*1.0D-5) / (2.0D0 * CPI) 
   end subroutine initialise_constants
   
end module constants



! Initialise constants, read in tables of physical data that are unchanged
! during the course of one run      
subroutine setsup
  use real_kind
  use constants
  use ltm2ubv
  
  implicit none
  
  ! Initialise some constants that cannot be computed at compile time on some
  ! compilers
  call initialise_constants
  
  ! Read nuclear reaction rates and neutrino loss parameters
  ! Used to be read from the same file as the opacity data, which is stupid
  !  since this is independent of metalicity
  call load_reaction_neutrino_rates(42)
  
  ! Read Bol. Corr, U-B, B-V table. 
  call load_colour_conversion_table(21)
  
  ! Read nuclear reaction (QRT) and neutrino (QNT) Q values, in MeV; constants 
  ! for electron screening (CZA, CZB, CZC, CZD, VZ); atomic parameters (CBN, KZN),
  ! with masses (CAN) consistent with Q-values; ionization potentials (CHI) and 
  ! statistical weights (COM); molecular hydrogen parameters (CH2)
  call load_atomic_data(26)
  
end subroutine setsup


