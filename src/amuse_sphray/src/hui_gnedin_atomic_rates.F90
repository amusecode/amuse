!> \file hui_gnedin_atomic_rates.F90

!> \brief Module that stores the rate fits from   
!! http://adsabs.harvard.edu/abs/1997MNRAS.292...27H
!!
!! Ferland3_1e9 are fits to the data from Ferland 1992 
!! and are accurate to 2% from 3K - 1e9K
!!
!! BurgessSeaton5e3_5e5 are fits to the data from Burgess
!! and Seaton 1960 and are accurate to 10% from 5e3K - 5e5K
!!
!! Lotz1e4_1e9 are fits to the data from Lotz 1967 
!! and are accurate to 3% from 1e4K to 1e9K
!!
!! AP3e4_1e6 are fits to the data from Aldrovandi and Pequignot 1973
!! and are accurate to 5% from 3e4K to 1e6K
!!
!! Black5e3_5e5 are from Black 1981 with the correction from Cen 1992
!! and are accurate to 10% from 5e3 to 5e5
!<
module hui_gnedin_atomic_rates_mod
use myf03_mod
implicit none

! contains
! HII   recomb case A, HII   recomb cool case A
! HeII  recomb case A, HeII  recomb cool case A
! HeIII recomb case A, HeIII recomb cool case A
! HII   recomb case B, HII   recomb cool case B
! HeII  recomb case B, HeII  recomb cool case B
! HeIII recomb case B, HeIII recomb cool case B
! He dielectronic recomb, He dielectronic recomb cool 

! HI   col ion, HI   col ion cool
! HeI  col ion, HeI  col ion cool
! HeII col ion, HeII col ion cool

! HI and HII col ext cool

  real(r8b), parameter, private :: h_erg_s = 6.626068d-27      !< Plank [ergs s]
  real(r8b), parameter, private :: h_eV_s = 4.1356668d-15      !< Plank [eV s]
  real(r8b), parameter, private :: k_erg_K = 1.3806503d-16     !< Boltzman [erg/K]
  real(r8b), parameter, private :: k_eV_K = 8.6173423d-5       !< Boltzman [eV/K]
  real(r8b), parameter, private :: nu_H = 13.6d0 / h_eV_s      !< HI   thresh [Hz]
  real(r8b), parameter, private :: nu_HeI = 24.587d0 / h_eV_s  !< HeI  thresh [Hz]
  real(r8b), parameter, private :: nu_HeII = 54.416d0 / h_eV_s !< HeII thresh [Hz] 
  real(r8b), parameter, private :: T_HI = 1.57807d5   !< HI   ion thresh in K
  real(r8b), parameter, private :: T_HeI = 2.85335d5  !< HeI  ion thresh in K
  real(r8b), parameter, private :: T_HeII = 6.31515d5 !< HeII ion thresh in K
 
contains

!--------------------------------------------------
!> HII recomb rate case A [cm^3/s] - Ferland3_1e9
  function Hui_HII_recombA(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (cm^3/s)
    real(r8b) :: lambda

    lambda = 2.d0 * T_HI / T
    rate = 1.269d-13 * lambda**1.503d0 / &
         ( 1.d0 + (lambda/0.522d0)**0.470d0 )**1.923d0

  end function Hui_HII_recombA

!-------------------------------------------------------------
!> HII recomb cooling rate case A [erg cm^3/s] - Ferland3_1e9
  function Hui_HII_rec_coolA(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (erg cm^3 / s)
    real(r8b) :: lambda

    lambda = 2.d0 * T_HI / T
    rate = 1.778d-29 * T * lambda**1.965d0 / &
         ( 1.d0 + (lambda/0.541d0)**0.502d0 )**2.697d0
  end function Hui_HII_rec_coolA

!-----------------------------------------------------------
!> HeII recomb rate case A [cm^3/s] - BurgessSeaton5e3_5e5
  function Hui_HeII_recombA(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (cm^3/s)
    real(r8b) :: lambda

    lambda = 2.d0 * T_HeI / T
    rate = 3.0d-14 * lambda**0.654d0 
  end function Hui_HeII_recombA

!----------------------------------------------------------------------
!> HeII recomb cooling rate case A [erg cm^3/s] - BurgessSeaton5e3_5e5
  function Hui_HeII_rec_coolA(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (erg cm^3 / s)
 
    rate = k_erg_K * T * Hui_HeII_recombA(T) 
  end function Hui_HeII_rec_coolA

!---------------------------------------------------
!> HeIII recomb rate case A [cm^3/s] - Ferland3_1e9
  function Hui_HeIII_recombA(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (cm^3/s)
    real(r8b) :: lambda

    lambda = 2.d0 * T_HeII / T
    rate = 2.0d0 * 1.269d-13 * lambda**1.503d0 / & 
         ( 1.d0 + (lambda/0.522d0)**0.470d0 )**1.923d0
  end function Hui_HeIII_recombA

!---------------------------------------------------------------
!> HeIII recomb cooling rate case A [erg cm^3/s] - Ferland3_1e9
  function Hui_HeIII_rec_coolA(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (erg cm^3 / s)
    real(r8b) :: lambda

    lambda = 2.d0 * T_HeII / T
    rate = 8.0d0 * 1.778d-29 * T * lambda**1.965d0 / & 
         ( 1.d0 + (lambda/0.541d0)**0.502d0 )**2.697d0
  end function Hui_HeIII_rec_coolA

!--------------------------------------------------
!> HII recomb rate case B [cm^3/s] - Ferland3_1e9
  function Hui_HII_recombB(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (cm^3/s)
    real(r8b) :: lambda

    lambda = 2.d0 * T_HI / T
    rate = 2.753d-14 * lambda**1.500d0 / &
         ( 1.d0 + (lambda/2.740d0)**0.407d0 )**2.242d0

  end function Hui_HII_recombB

!-------------------------------------------------------------
!> HII recomb cooling rate case B [erg cm^3/s] - Ferland3_1e9
  function Hui_HII_rec_coolB(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (erg cm^3 / s)
    real(r8b) :: lambda

    lambda = 2.d0 * T_HI / T
    rate = 3.435d-30 * T * lambda**1.970d0 / &
         ( 1.d0 + (lambda/2.250d0)**0.376d0 )**3.720d0
  end function Hui_HII_rec_coolB

!---------------------------------------------------------------
!> HeII recomb rate case B [cm^3 s^-1] - BurgessSeaton5e3_5e5
  function Hui_HeII_recombB(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (cm^3/s)
    real(r8b) :: lambda

    lambda = 2.d0 * T_HeI / T
    rate = 1.26d-14 * lambda**0.750d0 
  end function Hui_HeII_recombB

!---------------------------------------------------------------------
!> HeII recomb cooling rate case B [erg cm^3/s] - BurgessSeaton5e3_5e5
  function Hui_HeII_rec_coolB(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (erg cm^3 / s)
  
    rate = k_erg_K * T * Hui_HeII_recombB(T) 
  end function Hui_HeII_rec_coolB

!---------------------------------------------------
!> HeIII recomb rate case B [cm^3/s] - Ferland3_1e9
  function Hui_HeIII_recombB(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (cm^3/s)
    real(r8b) :: lambda
 
    lambda = 2.d0 * T_HeII / T
    rate = 2.0d0 * 2.753d-14 * lambda**1.500d0 / & 
         ( 1.d0 + (lambda/2.740d0)**0.407d0 )**2.242d0
  end function Hui_HeIII_recombB

!---------------------------------------------------------------
!> HeIII recomb cooling rate case B [erg cm^3/s] - Ferland3_1e9
  function Hui_HeIII_rec_coolB(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (cm^3/s)
    real(r8b) :: lambda

    lambda = 2.d0 * T_HeII / T
    rate = 8.0d0 * 3.435d-30 * T * lambda**1.970d0 / & 
        ( 1.d0 + (lambda/2.250d0)**0.376d0 )**3.720d0
  end function Hui_HeIII_rec_coolB


!---------------------------------------------------------------
!> HI collisional ionization rate [cm^3/s] - Lotz1e4_1e9
  function Hui_HI_col_ion(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (cm^3/s)
    real(r8b) :: lambda

    lambda = 2.d0 * T_HI / T
    rate = 21.11d0 * T**(-1.5) * exp(-lambda*0.5) *  & 
        lambda**(-1.089) / ( 1.d0 + (lambda/0.354d0)**0.874d0 )**1.101d0
  end function Hui_HI_col_ion

!---------------------------------------------------------------
!> HI collisional ionization cooloing rate [erg cm^3/s] - Lotz1e4_1e9
  function Hui_HI_col_ion_cool(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (erg cm^3 / s)

    rate = k_erg_K * T_HI * Hui_HI_col_ion(T) 
  end function Hui_HI_col_ion_cool


!---------------------------------------------------------------
!> HeI collisional ionization rate [cm^3/s] - Lotz1e4_1e9
  function Hui_HeI_col_ion(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (cm^3/s)
    real(r8b) :: lambda

    lambda = 2.d0 * T_HeI / T
    rate = 32.28d0 * T**(-1.5) * exp(-lambda*0.5) * lambda**(-1.146) / &
         ( 1.d0 + (lambda/0.416d0)**0.987d0 )**1.056d0
  end function Hui_HeI_col_ion

!---------------------------------------------------------------
!> HeI collisional ionization cooloing rate [erg cm^3/s] - Lotz1e4_1e9
  function Hui_HeI_col_ion_cool(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (erg cm^3 / s)
 
    rate = k_erg_K * T_HeI * Hui_HeI_col_ion(T) 
  end function Hui_HeI_col_ion_cool

!---------------------------------------------------------------
!> HeII collisional ionization rate [cm^3/s] - Lotz1e4_1e9
  function Hui_HeII_col_ion(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (cm^3/s)
    real(r8b) :: lambda

    lambda = 2.d0 * T_HeII / T
    rate = 19.65d0 * T**(-1.5) * exp(-lambda*0.5) * lambda**(-1.089) / &
         ( 1.d0 + (lambda/0.553d0)**0.735d0 )**1.275d0
  end function Hui_HeII_col_ion

!---------------------------------------------------------------
!> HeII collisional ionization cooloing rate [erg cm^3/s] - Lotz1e4_1e9
  function Hui_HeII_col_ion_cool(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (erg cm^3 / s)
  
    rate = k_erg_K * T_HeII * Hui_HeII_col_ion(T) 
  end function Hui_HeII_col_ion_cool


!----------------------------------------------
!> He dielectronic recombination [cm^3 s^-1] - AP3e4_1e6 
  function Hui_He_dielec_recomb(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (cm^3/s)
    real(r8b) :: lambda

    lambda = 2.0d0 * T_HeII  / T
    rate = 1.90e-3 * T**(-1.5e0) * exp(-0.75*lambda*0.5) * &
         ( 1.0e0 + 0.3e0 * exp(-0.15 * lambda * 0.5) )
  end function Hui_He_dielec_recomb
    
!---------------------------------------------------------
!>  He dielectronic recombination cooling [erg cm^3 s^-1]
  function Hui_He_dielec_recomb_cool(T) result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (erg cm^3 / s)

    rate = 0.75 * k_erg_K * T_HeII * Hui_He_dielec_recomb(T)            
  end function Hui_He_dielec_recomb_cool


!----------------------------------------------------------------------------
!> HI collisional excitation cooling [erg cm^3 s^-1] - Black5e3_5e5  
  function Hui_HI_col_ext_cool(T)  result(rate)
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (erg cm^3 / s)
    real(r8b) :: lambda

    lambda = 2.0d0 * T_HI  / T
    rate = 7.5e-19  * exp(-0.75*lambda*0.5) / (1.0 + sqrt(T*1.0e-5))
  end function Hui_HI_col_ext_cool


!----------------------------------------------------------------------------
!> HeII collisional excitation cooling [erg cm^3 s^-1] - Black5e3_5e5  
  function Hui_HeII_col_ext_cool(T) result(rate) 
    real(r8b) :: T     !< input temperature (K)
    real(r8b) :: rate  !< output rate (erg cm^3 / s)
    real(r8b) :: lambda

    lambda = 2.0d0 * T_HeII  / T
    rate = 5.54e-17 * (1.0d0/T)**(0.397) * &
         exp(-0.75*lambda*0.5) / (1.0 + sqrt(T*1.0e-5))
  end function Hui_HeII_col_ext_cool


end module hui_gnedin_atomic_rates_mod
