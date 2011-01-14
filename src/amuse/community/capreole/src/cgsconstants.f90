!>
!! \brief This module contains physical constants and conversion factors
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 
!!
!! \b Version: cgs units
!!

module cgsconstants

  use precision, only: dp
  use mathconstants, only: pi

  implicit none

  ! A collection of physical constants and conversion factors
  ! Units: cgs

  !> proton mass
  real(kind=dp), parameter :: m_p=1.672661e-24_dp
  !> speed of light
  real(kind=dp), parameter :: c=2.997925e+10_dp
  !> Planck constant
  real(kind=dp), parameter :: hplanck=6.6260755e-27_dp
  !> Stefan-Boltzmann constant
  real(kind=dp), parameter :: sigmasb=5.670e-5_dp
  !> Boltzmann constant
  real(kind=dp), parameter :: kb=1.381e-16_dp
  !> Gravitational constant
  real(kind=dp), parameter :: G_grav=6.6732e-8_dp

  !> ev2k   - conversion factor between evs and kelvins
  real(kind=dp),parameter  :: ev2k=1.0/8.617e-05
  !> ev2erg  - conversion factor between evs and ergs
  real(kind=dp),parameter  :: ev2erg=1.602e-12
  !> ev2j   - conversion factor between ergs and Joules
  real(kind=dp),parameter  :: erg2j=1e-7
  
  ! The following are scaled to frequency scaling

  !> Frequency scaling factor,
  !! this scaling parameter is independent of any main program scaling
  !! (see scaling.f90), and may only be used in the radiation physics 
  !! subroutines (currently switched off)
  real(kind=dp),parameter  :: sclfre=1.0e15
  !> conversion between evs and frequency
  real(kind=dp),parameter  :: ev2fr=0.241838e15!/sclfre      

  ! h/k, Planck/Boltzm
  ! Check this number, seems scaled
  !real(kind=dp),parameter  :: hoverk=47979.72484
  !> Planck constant scaled
  real(kind=dp),parameter  :: hscl=hplanck!*sclfre 
  !> tpic2  - 2*pi/c^2 times scaling factors needed for the integral cores
  real(kind=dp),parameter :: tpic2=2.0*pi/(c*c)
  !> two_pi_c2  - 2*pi/c^2 times scaling factors needed for the integral cores
  real(kind=dp),parameter :: two_pi_c2=2.0*pi/(c*c)!*sclfre**3

  !> Hydrogen recombination parameter (power law index)
  real(kind=dp),parameter :: albpow=-0.7_dp  !in the old code -0.79
  !> Hydrogen recombination parameter (value at 10^4 K)
  real(kind=dp),parameter :: bh00=2.59e-13_dp ! OTS value, alpha_B
  !> Helium0 recombination parameter (power law index)
  real(kind=dp), parameter :: alcpow=-0.672_dp
  !> Helium0 recombination parameter (value at 10^4 K)
  real(kind=dp), parameter :: bhe00=4.26e-13_dp !alpha_b+alpha_1
  !> Helium1 recombination parameter (value at 10^4 K)
  real(kind=dp), parameter :: bhe10=1.53e-12_dp !different in the book! 
  !here it was 5.91e-12 I replace with book value of 2*1.1e-12

  !> Hydrogen ionization energy (in eV)
  real(kind=dp), parameter :: eth0=13.598
  !> Hydrogen ionization energy (in erg)
  real(kind=dp),parameter :: hionen=eth0*ev2erg
  !> Hydrogen ionization energy expressed in K
  real(kind=dp),parameter :: temph0=eth0*ev2k
  !> Hydrogen collisional ionization parameter 1
  real(kind=dp),parameter :: xih0=1.0
  !> Hydrogen collisional ionization parameter 2
  real(kind=dp),parameter :: fh0=0.83
  !> Hydrogen collisional ionization parameter
  real(kind=dp),parameter :: colh0=1.3e-8*fh0*xih0/(eth0*eth0)


  !> Helium ionization energy (in eV)
  real(kind=dp), dimension(0:1), parameter :: ethe=(/24.587,54.416/)
  !> Helium ionization energy (in erg)
  real(kind=dp), dimension(0:1), parameter :: heionen=(/ethe(0)*ev2erg,ethe(1)*ev2erg/)
  !> Helium ionization energy expressed in K
  real(kind=dp), dimension(0:1), parameter :: temphe=(/ethe(0)*ev2k,ethe(1)*ev2k/)
  !> Helium collisional ionization parameter 1
  real(kind=dp),dimension(0:1),parameter :: xihe=(/2.0,1.0/)
  !> Helium collisional ionization parameter 2
  real(kind=dp),dimension(0:1),parameter :: fhe=(/0.63,1.30/)
  !> Helium collisional ionization parameter
  real(kind=dp),dimension(0:1),parameter :: colhe=(/1.3e-8*fhe(0)*xihe(0)/(ethe(0)*ethe(0)), &
  	1.3e-8*fhe(1)*xihe(1)/(ethe(1)*ethe(1))/)



end module cgsconstants


