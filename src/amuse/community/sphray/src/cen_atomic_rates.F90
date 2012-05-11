!> \file cen_atomic_rates.F90

!> \brief Module that stores the rate fits from  
!! http://adsabs.harvard.edu/abs/1992ApJS...78..341C
!<
module cen_atomic_rates_mod
use myf03_mod
implicit none

  real(r8b), parameter, private :: h_erg_s = 6.626068e-27 !< plank [erg s]
  real(r8b), parameter, private :: h_eV_s = 4.1356668e-15 !< plank [eV s]
  real(r8b), parameter, private :: k_erg_K = 1.3806503e-16 !< Boltzman [erg/K]
  real(r8b), parameter, private :: k_eV_K = 8.6173423e-5 !< Boltzman [eV/K]
  real(r8b), parameter, private :: nu_H = 13.6 / h_eV_s !< HI threshold [Hz] 
  real(r8b), parameter, private :: nu_HeI = 24.587 / h_eV_s !< HeI threshold [Hz]
  real(r8b), parameter, private :: nu_HeII = 54.416 / h_eV_s !< HeII threshold [Hz]
 
contains

! random stuff
!**************
! Osterbrok 1989
! HyRecA(5000)=6.82e-13 ; HyRecA(10000)=4.18e-13 ; HyRecA(20000)=2.51e-13
! HyRecB(5000)=4.54e-13 ; HyRecB(10000)=2.59e-13 ; HyRecB(20000)=2.52e-13
!

!#########################################!
! Heating and Cooling functions           !
!#########################################!


! collisional ionization cooling 
!================================
  
!----------------------------------------------------
!>  HI collisional ionization cooling [erg cm^3 s^-1]
    function Cen_HI_col_ion_cool(T)  
        real(r8b)  :: Cen_HI_col_ion_cool !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HI_col_ion_cool = 1.27e-21 * T**(1./2.) * &
                     (1 + (T/1.0e+5)**(1./2.))**(-1) * exp(-1.578091e+5/T)
    end function Cen_HI_col_ion_cool

!-----------------------------------------------------
!>  HeI collisional ionization cooling [erg cm^3 s^-1]
    function Cen_HeI_col_ion_cool(T)
        real(r8b)  :: Cen_HeI_col_ion_cool !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HeI_col_ion_cool = 9.38e-22 * T**(1./2.) * &
                     (1 + (T/1.0e+5)**(1./2.))**(-1) * exp(-2.853354e+5/T)
    end function Cen_HeI_col_ion_cool

!-----------------------------------------------------
!> HeII collisional ionization cooling [erg cm^3 s^-1]
    function Cen_HeII_col_ion_cool(T) 
        real(r8b)  :: Cen_HeII_col_ion_cool !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HeII_col_ion_cool = 4.95e-22 * T**(1./2.) * &
                     (1 + (T/1.0e+5)**(1./2.))**(-1) * exp(-6.31515e+5/T)
    end function Cen_HeII_col_ion_cool

!------------------------------------------------------
!> He(2^3S) collisional ionization cooling [cm^6 s^-1]
    function Cen_He23s_col_ion_cool(T)
        real(r8b)  :: Cen_He23s_col_ion_cool !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_He23s_col_ion_cool = 5.01e-27 * T**(-1.687e-1) * &
                     (1 + (T/1.0e+5)**(1./2.))**(-1) * exp(-5.5338e+4/T)
    end function Cen_He23s_col_ion_cool


! Recombination cooling 
!=========================

!--------------------------------------------
!> HII recombination cooling [erg cm^3 s^-1] 
    function Cen_HII_recomb_coolA(T)
        real(r8b)  :: Cen_HII_recomb_coolA !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HII_recomb_coolA = 8.70e-27 * T**(1./2.) * (T/1.0e+3)**(-0.2) * &
                     (1 + (T/1.0e+6)**(0.7))**(-1) 
    end function Cen_HII_recomb_coolA

!----------------------------------------------
!>   HeII recombination cooling [erg cm^3 s^-1]
    function Cen_HeII_recomb_coolA(T)
        real(r8b)  :: Cen_HeII_recomb_coolA !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HeII_recomb_coolA = 1.55e-26 * T**(0.3647e0)
    end function Cen_HeII_recomb_coolA

!-----------------------------------------------
!>   HeIII recombination cooling [erg cm^3 s^-1]
    function Cen_HeIII_recomb_coolA(T)
        real(r8b)  :: Cen_HeIII_recomb_coolA !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HeIII_recomb_coolA = 3.48e-26 * T**(1./2.) * (T/1.0e+3)**(-0.2) * &
                     (1 + (T/1.0e+6)**(0.7))**(-1) 
    end function Cen_HeIII_recomb_coolA

!---------------------------------------------------------
!>  He dielectronic recombination cooling [erg cm^3 s^-1]
    function Cen_He_dielec_recomb_cool(T)
        real(r8b)  :: Cen_He_dielec_recomb_cool !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_He_dielec_recomb_cool = 1.24e-13 * T**(-1.5e0) * &
                     exp(-4.7e5/T) * ( 1.0e0 + 0.3e0 * exp(-9.4e4/T) )
    end function Cen_He_dielec_recomb_cool


! Collisional excitation cooling 
!==================================

!-----------------------------------------------------------
!> HI collisional excitation cooling (all n) [erg cm^3 s^-1]
    function Cen_HI_col_ext_cool(T) 
        real(r8b)  :: Cen_HI_col_ext_cool !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HI_col_ext_cool = 7.50e-19 * &
                     (1 + (T/1.0e+5)**(1./2.))**(-1) * exp(-118348.e0/T)
    end function Cen_HI_col_ext_cool 

!-------------------------------------------------------------------------
!> HeI collisional excitation cooling (n = 2,3,4 triplets) [erg cm^6 s^-1]
    function Cen_HeI_col_ext_cool(T)
        real(r8b)  :: Cen_HeI_col_ext_cool !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HeI_col_ext_cool = 9.10e-27 * T**(-0.1687e0) * &
                     (1 + (T/1.0e+5)**(1./2.))**(-1) * exp(-13179.e0/T)
    end function Cen_HeI_col_ext_cool 

!---------------------------------------------------------------
!> HeII collisional excitation cooling (n = 2) [erg cm^3 s^-1]
    function Cen_HeII_col_ext_cool(T)
        real(r8b)  :: Cen_HeII_col_ext_cool !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HeII_col_ext_cool = 5.54e-17 * T**(-0.397e0) * &
                     (1 + (T/1.0e+5)**(1./2.))**(-1) * exp(-473638.e0/T)
    end function Cen_HeII_col_ext_cool 

!-----------------------------------------
!>  bremsstrahlung cooling [erg cm^3 s^-1]
    function Cen_Bremss_cool(T,nHII,nHeII,nHeIII,ne)
        real(r8b)  :: Cen_Bremss_cool !< rate
        real(r8b)  :: T !< temperature [K]
        real(r8b)  :: ne !< electron number density
        real(r8b)  :: nHII   !< HII number density
        real(r8b)  :: nHeII  !< HeII number density
        real(r8b) ::  nHeIII !< HeIII number density
        real(r8b)  :: gff !< gaunt factor
        gff = 1.5e0
        Cen_Bremss_cool = 1.42e-27 * sqrt(T) * gff * &
                          (nHII + nHeII + 4*nHeIII) * ne
    end function Cen_Bremss_cool

!-------------------------------------------------
!>  Haiman bremsstrahlung cooling [erg cm^3 s^-1]
    function Haiman_Bremss_cool_He(T,nHII,nHeII,nHeIII,ne)
        real(r8b)  :: Haiman_Bremss_cool_He !< rate
        real(r8b)  :: T !< temperature [K]
        real(r8b)  :: ne !< electron number density
        real(r8b)  :: nHII   !< HII number density
        real(r8b)  :: nHeII  !< HeII number density
        real(r8b) ::  nHeIII !< HeIII number density
        real(r8b)  :: gff !< gaunt factor
        gff = 1.10e0 + 0.34e0 * exp( -(5.50e0 - log10(T))**2/3.0e0 )
        Haiman_Bremss_cool_He = 1.42e-27 * T**(1./2.) * gff * &
                          (nHII + nHeII + 4*nHeIII) * ne
    end function Haiman_Bremss_cool_He

!-------------------------------------------------
!>  Haiman bremsstrahlung cooling [erg cm^3 s^-1]
    function Haiman_Bremss_cool_H(T,nHII,ne)
        real(r8b)  :: Haiman_Bremss_cool_H !< rate
        real(r8b)  :: T !< temperature [K]
        real(r8b)  :: ne !< electron number density
        real(r8b)  :: nHII   !< HII number density
        real(r8b)  :: gff !< gaunt factor
        gff = 1.10e0 + 0.34e0 * exp( -(5.50e0 - log10(T))**2/3.0e0 )
        Haiman_Bremss_cool_H = 1.42e-27 * T**(1./2.) * gff * &
                          (nHII) * ne
    end function Haiman_Bremss_cool_H

!---------------------------------------------------    
!> Haiman compton heating / cooling [erg cm^3 s^-1]
    function Haiman_Comp_heol(T,Tcmb,ne)
        real(r8b)  :: Haiman_Comp_heol !< rate
        real(r8b)  :: T     !< temperature [K]
        real(r8b)  :: Tcmb  !< background radiation temperature [K]
        real(r8b)  :: ne    !< electron number density
        Haiman_Comp_heol = 1.017e-37 * ne * Tcmb**4 * (T - Tcmb)  
    end function Haiman_Comp_heol


!#############################################!
! functions pertaining to the ionization state!
!#############################################!


! Collisional Ionization rates 
!==============================

!----------------------------------------
!> HI collisional ionization [cm^3 s^-1]   
    function Cen_HI_col_ion(T)
        real(r8b)  :: Cen_HI_col_ion !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HI_col_ion = 5.85e-11 * T**(1./2.) * &
                     (1.0 + (T/1.0e+5)**(1./2.))**(-1) * exp(-157809.1/T)
    end function Cen_HI_col_ion 

!----------------------------------------
!> HeI collisional ionization [cm^3 s^-1]  
    function Cen_HeI_col_ion(T)
        real(r8b)  :: Cen_HeI_col_ion !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HeI_col_ion = 2.38e-11 * T**(1./2.) * &
                     (1.0e+0 + (T/1.0e+5)**(1./2.))**(-1) * exp(-285335.4e0/T)
    end function Cen_HeI_col_ion 

!-----------------------------------------
!> HeII collisional ionization [cm^3 s^-1]  
    function Cen_HeII_col_ion(T)
        real(r8b)  :: Cen_HeII_col_ion !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HeII_col_ion = 5.68e-12 * T**(1./2.) * &
                     (1.0e+0 + (T/1.0e+5)**(1./2.))**(-1) * exp(-631515.e0/T)
    end function Cen_HeII_col_ion 


! Recombination rates, all from Cen 1992 except the one Spitzer
!===========================================================================

!---------------------------------------
!> HII recombination case A [cm^3 s^-1] 
    function Cen_HII_recombA(T)
        real(r8b)  :: Cen_HII_recombA !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HII_recombA = 8.40e-11 * T**(-1./2.) * (T/1.0e+3)**(-0.2) * &
                     (1 + (T/1.0e+6)**(0.7))**(-1)      
    end function Cen_HII_recombA


!---------------------------------------------------------------------
!> HII recombination case B (Power Law from Spitzer Table) [cm^3 s^-1] 
    function Spitzer_HII_recombB(T)
        real(r8b)  :: Spitzer_HII_recombB !< rate
        real(r8b)  :: T !< temperature [K]
        Spitzer_HII_recombB = 10**(-9.66286) / T**(0.739505)
    end function Spitzer_HII_recombB

!----------------------------------------
!> HeII recombination case A [cm^3 s^-1] 
    function Cen_HeII_recombA(T)
        real(r8b)  :: Cen_HeII_recombA !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_HeII_recombA = 1.50e-10 * T**(-0.6353e0)
    end function Cen_HeII_recombA

!-----------------------------------------
!> HeIII recombination case A [cm^3 s^-1]  
    function Cen_HeIII_recombA(T)
      real(r8b)  :: Cen_HeIII_recombA !< rate
      real(r8b)  :: T !< temperature [K]
      Cen_HeIII_recombA =  3.36e-10 * T**(-1./2.) * (T/1.0e+3)**(-0.2) * &
           ( 1 + (T/1.0e+6)**(0.7) )**(-1) 
    end function Cen_HeIII_recombA
    
!----------------------------------------------
!> He dielectronic recombination [cm^3 s^-1] 
    function Cen_He_dielec_recomb(T)
        real(r8b)  :: Cen_He_dielec_recomb !< rate
        real(r8b)  :: T !< temperature [K]
        Cen_He_dielec_recomb =  1.90e-3 * T**(-1.5e0) * &
                   exp(-470000.e0/T) * ( 1.0e0 + 0.3e0 * exp(-94000.e0/T) )
    end function Cen_He_dielec_recomb


! photoionization cross-sections from Verner
!===========================================================

  !> HI photo ionization x-section (Verner) [cm^2]
  !-------------------------------------------------------------------------
  function Verner_HI_photo_cs( Ry ) result( sigma )
    real(r8b), intent(in) :: Ry  !< energy [Rydbergs]
    real(r8b) :: sigma           !< cross section [cm^2]
    real(r8b), parameter :: Eth = 13.6d0
    real(r8b), parameter :: Emax = 5.0d4
    real(r8b), parameter :: E0 = 4.298d-1
    real(r8b), parameter :: sig0 = 5.475d4
    real(r8b), parameter :: ya = 3.288d1
    real(r8b), parameter :: P = 2.963d0

    real(r8b) :: eV
    real(r8b) :: x
    real(r8b) :: y

    eV = Ry * 13.6d0
    x = eV / E0
    y = x
  
    sigma = sig0 * (x-1)**2 * y**(0.5d0 * P - 5.5d0) * (1 + sqrt(y/ya))**(-P)
    sigma = sigma * 1.0d-18

  end function Verner_HI_photo_cs


! photoionization cross-sections from Osterbrok 1989
!===========================================================


!-----------------------------------------
!> HI photo ionization (Osterbrok) [cm^2]
    function Osterbrok_HI_photo_cs(freq) result(sigma)    
        real(r8b) :: freq  !< frequency [Hz]
        real(r8b) :: sigma !< cross section 
        if (freq < nu_H) then
           sigma = 0.0e0
           return
        end if
        sigma = 6.3e-18 * ( freq / nu_H ) ** (-3)
    end function Osterbrok_HI_photo_cs

!------------------------------------------
!> HeI photo ionization (Osterbrok)  [cm^2]
    function Osterbrok_HeI_photo_cs(freq) result(sigma)    
        real(r8b) :: freq  !< frequency [Hz]
        real(r8b) :: sigma !< cross section
        real(r8b) :: scaled_freq
        if (freq < nu_HeI) then
           sigma = 0.0e0
           return
        end if
        scaled_freq = freq / nu_HeI
        sigma = 7.2e-18 * & 
                ( 1.66e0 * ( scaled_freq )**(-2.05e0) + &
                  0.66e0 * ( scaled_freq )**(-3.05e0) )
    end function Osterbrok_HeI_photo_cs

!------------------------------------------
!> HeII photo ionization (Osterbrok)  [cm^2]
    function Osterbrok_HeII_photo_cs(freq) result(sigma)    
        real(r8b) :: freq  !< frequency
        real(r8b) :: sigma !< cross section
        if (freq < nu_HeII) then
           sigma = 0.0e0
           return
        end if
        sigma = 1.58e-18 * ( freq / nu_HeII ) ** (-3)
    end function Osterbrok_HeII_photo_cs


end module cen_atomic_rates_mod
