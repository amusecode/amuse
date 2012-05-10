!> \file bdf.F90

!> \brief the module for the photo backward difference solver.
!! this file contains the subroutines for updating particles
!! from photo ionizations, collisional ionizations, recombinations,
!! photo heating and atomic cooling.
!! the file bdf_np describes the routines to do the updates
!! without photo ionizations/heating.  
!<
module bdf_mod
use ionpar_mod
use cen_atomic_rates_mod
use atomic_rates_mod
use physical_constants_mod
use global_mod, only: rtable, isoT_k
implicit none

  integer(i8b), parameter :: MaxSteps = 50000000  !< maximum number of steps
  real(r8b), parameter :: ION_FAC = 0.01          !< ionization time factor
  real(r8b), parameter :: COOL_FAC = 0.01         !< cooling time factor
  real(r8b), parameter :: TIME_TOL = 0.001        !< fraction of time step to leave undone
 
  real(r8b), parameter :: zero = 0.0d0
  real(r8b), parameter :: one = 1.0d0
  real(r8b), parameter :: two = 2.0d0
  real(r8b), parameter :: three = 3.0d0

contains

!> this routine takes in the initial values in ipar 
!! {xHII,xHeII,xHeIII,T,pdep} and evolves the ionization particle
!! for a timestep ip%dt.  The final values are in the variables
!! ip%xHII, ip%xHeII, ip%xHeIII, ip%T, ip%pdeps
!========================================================================
subroutine bdfint(ip,scalls,photo,caseA,He,isoT,fixT)

  type(ionpart_type), intent(inout) :: ip !< ionization particle
  integer(i8b), intent(out) :: scalls !< number of runs through the solver
  logical, intent(in) :: photo   !< consider photo ionization / heating ?
  logical, intent(in) :: caseA(2)!< H/He case A recombinations ?
  logical, intent(in) :: He      !< consider Helium ?
  logical, intent(in) :: isoT    !< evolving temperature ? 
  logical, intent(in) :: fixT    !< fixed seperate temperature for each particle?

  type(atomic_rates_type) :: k
  real(r8b) :: t_tot
  real(r8b) :: dt2
  real(r8b) :: dt_i
  real(r8b) :: t_min
  integer(i8b) :: i
 
  ! check logicals
  !-------------------------------------------------------------------
  if (isoT .and. fixT) stop "isoT and fixT cant both be true"
  if (ip%dt_code <= zero) return


  ! initialize call counters and time variables
  !-------------------------------------------------------------------
  scalls = 0   
  t_tot = zero
  dt2 = ip%dt_s / two


  ! if constant temperature run, set constant atomic rates
  !-------------------------------------------------------------------
  if (isoT) then
     k = isoT_k
  end if

  if (fixT) then
     call get_atomic_rates(ip%T,rtable,k)
  end if


  ! do first round of calculations (only have xHI,xHeI,xHeII,nH,nHe)
  !-------------------------------------------------------------------
  ip%nHI  = ip%nH * ip%xHI
  ip%nHII = ip%nH * ip%xHII
  ip%ne = ip%nHII + ip%NeBckgnd

  if (He) then
     ip%nHeI   = ip%nHe * ip%xHeI
     ip%nHeII  = ip%nHe * ip%xHeII     
     ip%nHeIII = ip%nHe * ip%xHeIII
     ip%ne = ip%nHII + ip%nHeII + two * ip%nHeIII + ip%NeBckgnd
  end if

  if (photo) then
     call set_gammas(ip,He)
  else
     ip%gammaHI = zero
     ip%gammaHeI = zero
     ip%gammaHeII = zero
  end if


  ! begin sub-cycling steps
  !-------------------------------------------------------------------
  do i = 1,MaxSteps

     ip%iter = i


     ! calculate the cooling function
     !-------------------------------------------------------------------
     if (.not. isoT .and. .not. fixT) then
        ip%u = 1.5 * (ip%nH + ip%nHe + ip%ne) * k_erg_K * ip%T 
        call get_atomic_rates(ip%T,rtable,k)
        call set_cooling_func(ip,k,photo,caseA,He)
        if (photo) then
           ip%dudt = ip%COOLp
        else
           ip%dudt = ip%COOL
        end if
        ip%dTdt = two * ip%dudt / (three * (ip%nH + ip%nHe + ip%ne) * k_erg_K)
        if (abs(ip%dudt) > 0.0) then
           ip%tcool = abs(ip%u / ip%dudt)
        else
           ip%tcool = dt2
        endif
     end if


     ! calculate the time rate of change of the electron density
     !-------------------------------------------------------------------
     call set_ionization_func(ip,k,photo,caseA,He)
 
     call set_dnedt(ip,photo,caseA,He)
     if (abs(ip%dnedt) > 0.0) then
        ip%tion = abs(ip%ne / ip%dnedt)
     else
        ip%tion = dt2
     endif


     ! calculate the time step
     !-------------------------------------------------------------------
     if (isoT .or. fixT) then
        dt_i = ION_FAC * ip%tion
     else
        dt_i = min(ION_FAC * ip%tion, COOL_FAC * ip%tcool)
     end if
     
     if (dt_i > ip%dt_s - t_tot) dt_i = ip%dt_s - t_tot
     if (dt_i > dt2) dt_i = dt2


     ! take the time step
     !-------------------------------------------------------------------
     call take_bdf_ion_step(ip,k,dt_i,He)
     
     ip%strtag = "just after bdf_ion_step"
     call check_x(ip)


     if (.not. isoT .and. .not. fixT) then
        ip%T = ip%T + ip%dTdt * dt_i
        if (ip%T > 1.0d9) ip%T = 1.0d9
        if (ip%T < 1.0d0) ip%T = 1.0d0
     end if


     ! track photons, recombinations and update total time
     !-------------------------------------------------------------------

     !! photo ionizations

     if (photo) then
        ip%pdeps = ip%pdeps + ip%gammaHI * ip%HIcnt * dt_i
        if (He) then
           ip%pdeps = ip%pdeps + ip%gammaHeI  * ip%HeIcnt  * dt_i
           ip%pdeps = ip%pdeps + ip%gammaHeII * ip%HeIIcnt * dt_i
        end if
     end if


     !! time step

     t_tot = t_tot + dt_i
     t_min = huge(one)
     if (t_tot < t_min) t_min = t_tot
     if (abs(ip%dt_s - t_min) < TIME_TOL * ip%dt_s) exit
    
     if (i == MaxSteps) then
        write(*,*) " ************* BDF ERROR ************* "
        write(*,*) "finished ", MaxSteps, " little steps before "
        write(*,*) "completing large timestep.  "
        write(*,*) "t/dt = ", t_tot/dt_i, "rayn = ", ip%rayn
        write(*,*) "ipar = "
        call ionpar2screen(ip)
        write(*,*) " *************  bdf.f90  *************"
        stop
     end if

     scalls = scalls + 1
     if(photo) call set_gammas(ip,He)

  end do


  
end subroutine bdfint
 

!> this routine evolves the ionization fractions and possibly the 
!! temperature of the ionization particle forward in time by dt
!================================================================
subroutine take_bdf_ion_step(ip,k,dt,He)
use physical_constants_mod, only: k_erg_K

  type(ionpart_type), intent(inout) :: ip   !< ionization particle  
  type(atomic_rates_type), intent(in) :: k  !< rates
  real(r8b), intent(inout) :: dt                 !< time step
  logical, intent(in) :: He

  real(r8b) :: CC, DD
  real(r8b) :: nHIp, nHIIp, nep
  real(r8b) :: nHeIp, nHeIIp, nHeIIIp  
  real(r8b) :: mfac

  
  ! HI
  CC = ip%nHII * ip%ne * k%HIIrcB
  DD = ip%ne * k%HIci 
  DD = DD + ip%gammaHI
  nHIp = ( CC * dt + ip%nHI ) / (1.0d0 + DD * dt)

  ! HII
  CC = nHIp * ip%ne * k%HIci
  CC = CC + nHIp * ip%gammaHI
  DD = ip%ne * k%HIIrccB
  nHIIp = ( CC * dt + ip%nHII ) / (1.0d0 + DD * dt)


  if (He) then
    
     ! HeI
     CC = ip%nHeII * ip%ne * k%HeIIrcB
     DD = ip%ne * k%HeIci 
     DD = DD + ip%gammaHeI 
     nHeIp = ( CC * dt + ip%nHeI ) / (1.0d0 + DD * dt)
     
     ! HeII
     CC = (nHeIp * k%HeIci + ip%nHeIII * k%HeIIIrcB) * ip%ne
     DD = ip%ne * (k%HeIIrcB + k%HeIIci) 
     CC = CC + nHeIp * ip%gammaHeI
     DD = DD + ip%gammaHeII
     nHeIIp = ( CC * dt + ip%nHeII ) / (1.0d0 + DD * dt)
     
     ! HeIII
     CC = nHeIIp * ip%ne * k%HeIIci
     CC = CC + nHeIIp * ip%gammaHeII
     DD = ip%ne * k%HeIIIrcB
     nHeIIIp = ( CC * dt + ip%nHeIII ) / (1.0d0 + DD * dt)
     
  end if

  ! electrons
  CC = ip%gammaHI * nHIp    
  DD = k%HIIrcB   * nHIIp - k%HIci   * nHIp    

  if (He) then
     CC = CC &
          + ip%gammaHeI  * nHeIp   &
          + ip%gammaHeII * nHeIIp
     DD = DD &
          + k%HeIIrcB  * nHeIIp  - k%HeIci  * nHeIp   &
          + k%HeIIIrcB * nHeIIIp - k%HeIIci * nHeIIp
  end if

  nep   = ( CC * dt + ip%ne) / (1.0d0 + DD * dt)


  ! mass corrections

  mfac = ip%nH / (nHIp + nHIIp)
  ip%nHI  = nHIp  * mfac
  ip%nHII = nHIIp * mfac
  ip%xHI  = ip%nHI  / ip%nH
  ip%xHII = ip%nHII / ip%nH
  ip%ne = ip%nHII

  if (He) then
     mfac = ip%nHe / (nHeIp + nHeIIp + nHeIIIp)
     ip%nHeI   = nHeIp   * mfac
     ip%nHeII  = nHeIIp  * mfac
     ip%nHeIII = nHeIIIp * mfac
     ip%xHeI   = ip%nHeI   / ip%nHe
     ip%xHeII  = ip%nHeII  / ip%nHe
     ip%xHeIII = ip%nHeIII / ip%nHe
     ip%ne = ip%ne + ip%nHeII + 2.0d0 * ip%nHeIII
  end if

  ip%ne = ip%ne + ip%NeBckgnd

 
end subroutine take_bdf_ion_step

end module bdf_mod



