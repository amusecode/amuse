!> \file hummer_atomic_rates.F90

!> \brief Module that stores the rate table from  
!! http://adsabs.harvard.edu/abs/1994MNRAS.268..109H
!<
module hummer_atomic_rates_mod
use myf03_mod
implicit none
private

public :: Hum_HII_recombA, Hum_HII_recombB
public :: Hum_HII_recomb_coolA, Hum_HII_recomb_coolB

real(r8b), parameter :: k_erg_K = 1.3806d-16    !< Boltzman's (erg/K)

real(r8b), parameter :: logTmin = 1.0d0
real(r8b), parameter :: logTmax = 7.0d0

real(r8b), parameter :: Tmin = 10.0d0**(int(logTmin))
real(r8b), parameter :: Tmax = 10.0d0**(int(logTmax))

real(r8b), parameter :: dlogT = 0.2d0
integer(i8b), parameter :: Ntable = 31

real(r8b), parameter :: LogTemp_table(Ntable) =              &
  (/                                                         &
  1.0d0, 1.2d0, 1.4d0, 1.6d0, 1.8d0,      &
  2.0d0, 2.2d0, 2.4d0, 2.6d0, 2.8d0,      &
  3.0d0, 3.2d0, 3.4d0, 3.6d0, 3.8d0,      &
  4.0d0, 4.2d0, 4.4d0, 4.6d0, 4.8d0,      &
  5.0d0, 5.2d0, 5.4d0, 5.6d0, 5.8d0,      &
  6.0d0, 6.2d0, 6.4d0, 6.6d0, 6.8d0,      &
  7.0d0 /)

real(r8b), parameter :: HIIrc1_table(Ntable) =               &
  (/                                                         &
  1.646d-11, 1.646d-11, 1.646d-11, 1.646d-11, 1.646d-11,     &
  1.646d-11, 1.645d-11, 1.645d-11, 1.644d-11, 1.642d-11,     &
  1.640d-11, 1.636d-11, 1.629d-11, 1.620d-11, 1.605d-11,     &
  1.582d-11, 1.548d-11, 1.499d-11, 1.431d-11, 1.341d-11,     &
  1.227d-11, 1.093d-11, 9.454d-12, 7.920d-12, 6.427d-12,     &
  5.058d-12, 3.866d-12, 2.877d-12, 2.089d-12, 1.485d-12,     &
  1.036d-12 /)


real(r8b), parameter :: HIIrcB_table(Ntable) =                &  
  (/                                                          &
  9.283d-11,  8.823d-11,  8.361d-11,  7.898d-11,  7.435d-11,  &
  6.973d-11,  6.512d-11,  6.054d-11,  5.599d-11,  5.147d-11,  &
  4.700d-11,  4.258d-11,  3.823d-11,  3.397d-11,  2.983d-11,  &
  2.584d-11,  2.204d-11,  1.847d-11,  1.520d-11,  1.226d-11,  &
  9.696d-12,  7.514d-12,  5.710d-12,  4.257d-12,  3.117d-12,  &
  2.244d-12,  1.590d-12,  1.110d-12,  7.642d-13,  5.199d-13,  &
  3.498d-13  /)      

real(r8b), parameter :: HIIrcA_table(Ntable) = HIIrc1_table + HIIrcB_table


real(r8b), parameter :: HIIrcc1_table(Ntable) =               &
  (/                                                          &
  1.646d-11,  1.646d-11,  1.646d-11,  1.646d-11,  1.646d-11,  &
  1.645d-11,  1.644d-11,  1.643d-11,  1.641d-11,  1.638d-11,  &
  1.633d-11,  1.625d-11,  1.613d-11,  1.594d-11,  1.565d-11,  &
  1.522d-11,  1.460d-11,  1.374d-11,  1.260d-11,  1.119d-11,  &
  9.571d-12,  7.844d-12,  6.146d-12,  4.601d-12,  3.295d-12,  &
  2.262d-12,  1.494d-12,  9.520d-13,  5.878d-13,  3.528d-13,  &
  2.066d-13   /)

real(r8b), parameter :: HIIrccB_table(Ntable) =               &
  (/                                                          &
  8.287d-11,  7.821d-11,  7.356d-11,  6.892d-11,  6.430d-11,  &
  5.971d-11,  5.515d-11,  5.062d-11,  4.614d-11,  4.170d-11,  &
  3.734d-11,  3.306d-11,  2.888d-11,  2.484d-11,  2.098d-11,  &
  1.736d-11,  1.402d-11,  1.103d-11,  8.442d-12,  6.279d-12,  &
  4.539d-12,  3.192d-12,  2.185d-12,  1.458d-12,  9.484d-13,  &
  6.023d-13,  3.738d-13,  2.268d-13,  1.348d-13,  7.859d-14,  &
  4.499d-14 /)

real(r8b), parameter :: HIIrccA_table(Ntable) = HIIrcc1_table + HIIrccB_table

contains

!---------------------------------------
!> HII recombination case A [cm^3 s^-1] 
    function Hum_HII_recombA(T) result(rate)
        real(r8b)  :: T    !< temperature [K]
        real(r8b)  :: rate !< rate

        real(r8b) :: Tnum, Tremain, rdif
        integer(i8b) :: Tindx

        
        if (T < Tmin) stop "T < Tmin in hummer_atomic_rates_mod.f90"
    
        Tnum = ( log10(T) - logTmin ) / dlogT
        Tindx = ceiling(Tnum)
    
        if (Tindx == 0) then
           Tindx = 1
           Tremain = 0.0d0
        end if
        Tremain = Tnum - (Tindx-1)
    
        if(Tindx < 1)        stop "Tindx < 1 in hummer_atomic_rates_mod.f90 "
        if(Tindx > Ntable-1) stop "Tindx > Ntable-1 in hummer_atomic_rates_mod.f90 "

        rdif = (HIIrcA_table(Tindx+1) - HIIrcA_table(Tindx)) * Tremain
        rate = (HIIrcA_table(Tindx) + rdif) / sqrt(T)

    end function Hum_HII_recombA


!---------------------------------------
!> HII recombination case B [cm^3 s^-1] 
    function Hum_HII_recombB(T) result(rate)
        real(r8b)  :: T    !< temperature [K]
        real(r8b)  :: rate !< rate

        real(r8b) :: Tnum, Tremain, rdif
        integer(i8b) :: Tindx

        
        if (T < Tmin) stop "T < Tmin in hummer_atomic_rates_mod.f90"
    
        Tnum = ( log10(T) - logTmin ) / dlogT
        Tindx = ceiling(Tnum)
    
        if (Tindx == 0) then
           Tindx = 1
           Tremain = 0.0d0
        end if
        Tremain = Tnum - (Tindx-1)
    
        if(Tindx < 1)        stop "Tindx < 1 in hummer_atomic_rates_mod.f90 "
        if(Tindx > Ntable-1) stop "Tindx > Ntable-1 in hummer_atomic_rates_mod.f90 "

        rdif = (HIIrcB_table(Tindx+1) - HIIrcB_table(Tindx)) * Tremain
        rate = (HIIrcB_table(Tindx) + rdif) / sqrt(T)

    end function Hum_HII_recombB


!---------------------------------------
!> HII recombination cool case A [cm^3 s^-1] 
    function Hum_HII_recomb_coolA(T) result(rate)
        real(r8b)  :: T    !< temperature [K]
        real(r8b)  :: rate !< rate

        real(r8b) :: Tnum, Tremain, rdif
        integer(i8b) :: Tindx

        
        if (T < Tmin) stop "T < Tmin in hummer_atomic_rates_mod.f90"
    
        Tnum = ( log10(T) - logTmin ) / dlogT
        Tindx = ceiling(Tnum)
    
        if (Tindx == 0) then
           Tindx = 1
           Tremain = 0.0d0
        end if
        Tremain = Tnum - (Tindx-1)
    
        if(Tindx < 1)        stop "Tindx < 1 in hummer_atomic_rates_mod.f90 "
        if(Tindx > Ntable-1) stop "Tindx > Ntable-1 in hummer_atomic_rates_mod.f90 "

        rdif = (HIIrccA_table(Tindx+1) - HIIrccA_table(Tindx)) * Tremain
        rate = (HIIrccA_table(Tindx) + rdif) * sqrt(T) * k_erg_K

    end function Hum_HII_recomb_coolA


!---------------------------------------
!> HII recombination cool case B [cm^3 s^-1] 
    function Hum_HII_recomb_coolB(T) result(rate)
        real(r8b)  :: T    !< temperature [K]
        real(r8b)  :: rate !< rate

        real(r8b) :: Tnum, Tremain, rdif
        integer(i8b) :: Tindx

        
        if (T < Tmin) stop "T < Tmin in hummer_atomic_rates_mod.f90"
    
        Tnum = ( log10(T) - logTmin ) / dlogT
        Tindx = ceiling(Tnum)
    
        if (Tindx == 0) then
           Tindx = 1
           Tremain = 0.0d0
        end if
        Tremain = Tnum - (Tindx-1)
    
        if(Tindx < 1)        stop "Tindx < 1 in hummer_atomic_rates_mod.f90 "
        if(Tindx > Ntable-1) stop "Tindx > Ntable-1 in hummer_atomic_rates_mod.f90 "

        rdif = (HIIrccB_table(Tindx+1) - HIIrccB_table(Tindx)) * Tremain
        rate = (HIIrccB_table(Tindx) + rdif) * sqrt(T) * k_erg_K 

    end function Hum_HII_recomb_coolB

end module hummer_atomic_rates_mod
