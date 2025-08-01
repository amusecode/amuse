!> \file physical_constants.F90

!> \brief stores physical constants 
!!
!<

!> physical constants module
module physical_constants_mod
use myf03_mod
implicit none

! famous numbers
!-------------------

    real(r8b), parameter :: ex = 2.71828182846d0    !< base of the natural log
    real(r8b), parameter :: pi = 3.141592653589793238462643383279502884197d0 !< pi


! physical constants in cgs units (set to be compatible with Gadget)
!-------------------------------------------------------------------

    real(r8b), parameter :: h_erg_s = 6.6262d-27    !< Plank's const. (erg*s)
    real(r8b), parameter :: h_eV_s  = 4.1357d-15    !< Plank's const. (eV*s)
    real(r8b), parameter :: k_erg_K = 1.3806d-16    !< Boltzman's (erg/K)
    real(r8b), parameter :: k_eV_K  = 8.6173d-5     !< Boltzman's (eV/K)    
    real(r8b), parameter :: c = 2.9979d10           !< speed of light (cm/s)
    real(r8b), parameter :: G = 6.672d-8            !< grav const. (cm^3 g^-1 s^-2)

! usefull physical scales
!-------------------------

    real(r8b), parameter :: M_H  = 1.6737d-24       !< Hydrogen mass (g)
    real(r8b), parameter :: M_He = 6.6465d-24       !< Helium mass (g)  
    real(r8b), parameter :: M_p  = 1.6726d-24       !< proton mass (g) 
    real(r8b), parameter :: M_e  = 9.10953d-28      !< electron mass (g) 
    real(r8b), parameter :: M_solar = 1.989d33      !< solar mass (g)
    real(r8b), parameter :: L_solar = 3.826d33      !< solar luminosity (ergs/s) 
    real(r8b), parameter :: Rydberg = 1.36d1        !< ionization energy of HI (eV)
    real(r8b), parameter :: CMBtempNow = 2.728d0    !< CMB temperature now (K)
    
! convenient conversion factors
!------------------------------

    real(r8b), parameter :: Mpc2cm       = 3.085678d24         !< Mpc to cm
    real(r8b), parameter :: kpc2cm       = 3.085678d21         !< kpc to cm
    real(r8b), parameter :: pc2cm        = 3.085678d18         !< pc to cm
    real(r8b), parameter :: km2cm        = 1.0d5               !< km to cm

    real(r8b), parameter :: kpc2cm3      = kpc2cm**3           !< kpc^3 to cm^3
    real(r8b), parameter :: cm2kpc       = 3.24077649d-22      !< cm to kpc
    real(r8b), parameter :: cm2kpc3      = cm2kpc**3           !< cm^3 to kpc^3

    real(r8b), parameter :: Mpc2km       = 3.08568025d19       !< Mpc to km
    real(r8b), parameter :: km2Mpc       = 3.24077649d-20      !< km to Mpc
    real(r8b), parameter :: erg2eV       = 6.24150974d11       !< ergs to eV
    real(r8b), parameter :: eV2erg       = 1.60217646d-12      !< eV to ergs
    real(r8b), parameter :: erg2solar    = 1.0d0/L_solar       !< cgs 2 solar 
    real(r8b), parameter :: year2sec     = 3.1556926d7         !< years to seconds
    real(r8b), parameter :: Myr2sec      = 3.1556926d13        !< Myrs to seconds
    real(r8b), parameter :: Gyr2sec      = 3.1556926d16        !< Gyrs to seconds
    real(r8b), parameter :: sec2year     = 1.0d0 / year2sec    !< seconds to years
    real(r8b), parameter :: s2Myr        = 1.0d0 / Myr2sec     !< seconds to Myrs
    real(r8b), parameter :: g2solar      = 1.0d0 / M_solar     !< grams 2 solar masses
    real(r8b), parameter :: solarLMyrs   = L_solar * Myr2sec   !< cgs to ergs / Myrs
    real(r8b), parameter :: fourthirdspi = 4.188790204786391d0 !< four thirds * pi

! hydrogen and helium lines expressed in various ways
!---------------------------------------------------------------------

    real(r8b), parameter :: HI_th_eV   = 13.606d0  !< HI ionization threshold (eV)
    real(r8b), parameter :: HeI_th_eV  = 24.587d0  !< HeI ionization threshold (eV)
    real(r8b), parameter :: HeII_th_eV = 54.416d0  !< HeII ionization threshold (eV)
    
    real(r8b), parameter :: HI_th_erg   = HI_th_eV * eV2erg    !< HI ionization threshold (ergs)
    real(r8b), parameter :: HeI_th_erg  = HeI_th_eV * eV2erg   !< HeI ionization threshold (ergs)
    real(r8b), parameter :: HeII_th_erg = HeII_th_eV * eV2erg  !< HeII ionization threshold (ergs)

    real(r8b), parameter :: HI_th_Hz   = HI_th_eV   / h_eV_s  !< HI ionization threshold (Hz)
    real(r8b), parameter :: HeI_th_Hz  = HeI_th_eV  / h_eV_s  !< HeI ionization threshold (Hz)
    real(r8b), parameter :: HeII_th_Hz = HeII_th_eV / h_eV_s  !< HeII ionization threshold (Hz)
    
    real(r8b), parameter :: HI_th_cm   = h_eV_s * c / HI_th_eV    !< HI ionization threshold (cm)
    real(r8b), parameter :: HeI_th_cm  = h_eV_s * c / HeI_th_eV   !< HeI ionization threshold (cm)
    real(r8b), parameter :: HeII_th_cm = h_eV_s * c / HeII_th_eV  !< HeII ionization threshold (cm)

    real(r8b), parameter :: HI_th_nm   = HI_th_cm * 1.0d7    !< HI ionization threshold (nm)
    real(r8b), parameter :: HeI_th_nm  = HeI_th_cm * 1.0d7   !< HeI ionization threshold (nm)
    real(r8b), parameter :: HeII_th_nm = HeII_th_cm * 1.0d7  !< HeII ionization threshold (nm)

    real(r8b), parameter :: HI_th_A   = HI_th_cm * 1.0d8    !< HI ionization threshold (A)
    real(r8b), parameter :: HeI_th_A  = HeI_th_cm * 1.0d8   !< HeI ionization threshold (A)
    real(r8b), parameter :: HeII_th_A = HeII_th_cm * 1.0d8  !< HeII ionization threshold (A)

    real(r8b), parameter :: H_LyA_eV = HI_th_eV - HI_th_eV / 4.0d0  !< H Lyman Alpha (eV)
    real(r8b), parameter :: H_LyA_cm = h_eV_s * c / H_LyA_eV        !< H Lyman Alpha (cm)
    real(r8b), parameter :: H_LyA_nm = H_LyA_cm * 1.0d7             !< H Lyman Alpha (nm)
    real(r8b), parameter :: H_LyA_A  = H_LyA_cm * 1.0d8             !< H Lyman Alpha (A)
    real(r8b), parameter :: H_LyA_Hz = c / H_LyA_cm                 !< H Lyman Alpha (Hz)




end module physical_constants_mod
