!> \file euler.F90
 
!> \brief The module for the implicit euler solver.
!!
!<

module euler_mod
use ionpar_mod
use atomic_rates_mod
use physical_constants_mod
use global_mod, only: isoT_k, rtable

implicit none

  integer(i8b), parameter :: MaxSteps = 50000000 !< maximum number of steps
  real(r8b), parameter :: ION_FAC  = 1.0d-3      !< ionization time factor
  real(r8b), parameter :: COOL_FAC = 1.0d-3      !< cooling time factor
  real(r8b), parameter :: TIME_TOL = 1.0d-6      !< max frac. of step undone

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
subroutine eulerint(ip,scalls,photo,caseA,He,isoT,fixT)

  type(ionpart_type), intent(inout) :: ip !< ionization particle
  integer(i8b), intent(out) :: scalls !< number of runs through the solver
  logical, intent(in) :: photo   !< consider photo ionization / heating ?
  logical, intent(in) :: caseA(2)!< H/He case A recombinations ?
  logical, intent(in) :: He      !< consider Helium ?
  logical, intent(in) :: isoT    !< evolving temperature ? 
  logical, intent(in) :: fixT    !< fixed seperate temperature for each particle?

  real(r8b) :: xH(2)
  real(r8b) :: xHe(3)
  real(r8b) :: t_tot
  real(r8b) :: dt2
  real(r8b) :: dt_i
  real(r8b) :: t_min
  integer(i8b) :: i
  real(r8b) :: mfac
  type(atomic_rates_type) :: k

  ! initialize solver calls counter and time step
  !------------------------------------------------
  scalls = 0   
  t_tot  = zero
  dt2    = ip%dt_s / two

  ! if constant temperature run, set constant atomic rates
  !---------------------------------------------------------
  if (isoT .and. fixT) stop "isoT and fixT can't both be true"

  if (ip%dt_code.LE.0) return

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
  ip%ne   = ip%nHII + ip%NeBckgnd

  if (He) then
     ip%nHeI   = ip%nHe * ip%xHeI
     ip%nHeII  = ip%nHe * ip%xHeII     
     ip%nHeIII = ip%nHe * ip%xHeIII
     ip%ne = ip%nHII + ip%nHeII + two * ip%nHeIII + ip%NeBckgnd
  else
     ip%nHeI = zero
     ip%nHeII = zero
     ip%nHeIII = zero
  end if

  if (photo) then
     call set_gammas(ip,He)
  else
     ip%gammaHI   = zero
     ip%gammaHeI  = zero
     ip%gammaHeII = zero
  end if

!---------------------------------------------------------  
  
  do i = 1,MaxSteps

     ip%iter = i

     !----------------------------------------------------------  
     ! calculate the cooling function
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


     !----------------------------------------------------------  
     ! calculate the time rate of change of the electron density
     call set_ionization_func(ip,k,photo,caseA,He)
     call setDH(ip,photo,caseA(1))
     if (He) call setDHe(ip,photo,caseA(2))
     call set_dnedt(ip,photo,caseA,He)
     if (abs(ip%dnedt) > 0.0) then
        ip%tion = abs(ip%ne / ip%dnedt)
     else
        ip%tion = dt2
     endif


     !----------------------------------------------------------  
     ! calculate the time step
     if (isoT .or. fixT) then
        dt_i = ION_FAC * ip%tion
     else
        dt_i = min(ION_FAC * ip%tion, COOL_FAC * ip%tcool)
     end if
     
     if (dt_i > ip%dt_s - t_tot) dt_i = ip%dt_s - t_tot
     if (dt_i > dt2) dt_i = dt2

     !----------------------------------------------------------  
     ! take the time step

     xH(1) = ip%xHI
     xH(2) = ip%xHII
     xHe(1) = ip%xHeI
     xHe(2) = ip%xHeII
     xHe(3) = ip%xHeIII

     if (He) then
        call take_euler_ion_step(xH,ip%DH,dt_i,xHe,ip%DHe)
     else
        call take_euler_ion_step(xH,ip%DH,dt_i)
     end if

     ip%xHI  = xH(1)
     ip%xHII = xH(2)
     ip%xHeI   = xHe(1)
     ip%xHeII  = xHe(2)
     ip%xHeIII = xHe(3)

!     ip%strtag = "just after euler_ion_step"
!     call check_x(ip)

     !----------------------------------------------------------  
     ! do mass conservation

     mfac = one / (ip%xHI + ip%xHII)
     ip%xHI = ip%xHI * mfac
     ip%xHII = ip%xHII * mfac

     if (He) then
        mfac = one / (ip%xHeI + ip%xHeII + ip%xHeIII)
        ip%xHeI   = ip%xHeI   * mfac
        ip%xHeII  = ip%xHeII  * mfac
        ip%xHeIII = ip%xHeIII * mfac
     end if

     ip%strtag = "just after mass correction"
     call check_x(ip)


     ip%nHI  = ip%nH * ip%xHI
     ip%nHII = ip%nH * ip%xHII
     ip%ne   = ip%nHII + ip%NeBckgnd
 
     if (He) then
        ip%nHeI   = ip%nHe * ip%xHeI
        ip%nHeII  = ip%nHe * ip%xHeII
        ip%nHeIII = ip%nHe * ip%xHeIII
        ip%ne = ip%nHII + ip%nHeII + two * ip%nHeIII + ip%NeBckgnd
     end if


     if (.not. isoT .and. .not. fixT) then
        ip%T = ip%T + ip%dTdt * dt_i
        if (ip%T > 1.0d9) ip%T = 1.0d9
        if (ip%T < 1.0d0) ip%T = 1.0d0
!        write(*,*) "isoT,fixT", isoT,fixT
!        write(*,*) "no one is supposed to be here"
!        stop
!        ip%u = 1.5d0 * (ip%nH + ip%nHe + ip%ne) * k_erg_K * ip%T
!        call take_euler_temp_step(ip,dt_i)
     end if


     !----------------------------------------------------------  
     ! track photons, recombinations and update total time

     ! photo ionizations

     if (photo) then
        ip%pdeps = ip%pdeps + ip%gammaHI * ip%HIcnt * dt_i
        if (He) then
           ip%pdeps = ip%pdeps + ip%gammaHeI  * ip%HeIcnt  * dt_i
           ip%pdeps = ip%pdeps + ip%gammaHeII * ip%HeIIcnt * dt_i
        end if
     end if


     ! time step

     t_tot = t_tot + dt_i
     t_min = huge(one)
     if (t_tot < t_min) t_min = t_tot
     if (abs(ip%dt_s - t_min) < TIME_TOL * ip%dt_s) exit
    
     if (i == MaxSteps) then
        write(*,*) " ************* EULER ERROR ************* "
        write(*,*) "finished ", MaxSteps, " little steps before "
        write(*,*) "completing large timestep.  "
        write(*,*) "t/dt = ", t_tot/dt_i, "rayn = ", ip%rayn
        write(*,*) "ipar = "
        call ionpar2screen(ip)
        write(*,*) " *************  euler.f90  *************"
        stop
     end if

     scalls = scalls + 1
     if(photo) call set_gammas(ip,He)

  end do

  
end subroutine eulerint
 

!> this routine takes in the initial values in ipar 
!! {xHII,xHeII,xHeIII,T,pdep} and simply deposits recombination 
!! photons based on the optical depth.  The final values are in the variables
!! ip%xHII, ip%xHeII, ip%xHeIII, ip%T, ip%pdeps
!========================================================================
subroutine recombeulerint(ip,scalls)


  real(r8b), parameter :: TAU_TOL = 1.0e-4

  type(ionpart_type), intent(inout) :: ip !< ionization particle
  integer(i8b), intent(out) :: scalls !< number of runs through the solver
 
  scalls = 0   

! warning about the following:
  if (ip%dt_code.LE.0) return

!---------------------------------------------------------
! set initial number and ionization fractions

  ip%nHII = ip%nH * ip%xHII
  ip%nHI  = ip%nH - ip%nHII
  ip%xHI  = ip%nHI / ip%nH
  ip%ne = ip%nHII + ip%NeBckgnd
#ifdef incHe
  ip%nHeIII = ip%nHe * ip%xHeIII
  ip%nHeII  = ip%nHe * ip%xHeII
  ip%nHeI   = ip%nHe - ip%nHeIII - ip%nHeII
  ip%xHeI  = ip%nHeI / ip%nHe
  ip%ne = ip%nHII + ip%nHeII + 2.0 * ip%nHeIII + ip%NeBckgnd
#endif  

!---------------------------------------------------------
! simlple count of atoms that are absorbing photons

  ip%HIcnt = ip%Hcnt * ( ip%nHI / ip%nH )       ! N HI atoms
  ip%Allcnt = ip%HIcnt
#ifdef incHe
  ip%HeIcnt = ip%Hecnt * ( ip%nHeI / ip%nHe )   ! N HeI atoms
  ip%HeIIcnt = ip%Hecnt * ( ip%nHeII / ip%nHe ) ! N HeII atoms
  if (ip%penrg .GE. HeI_th_erg) then
     ip%Allcnt = ip%HIcnt + ip%HeIcnt
  else if (ip%penrg .GE. HeII_th_erg) then
     ip%Allcnt = ip%HIcnt + ip%HeIcnt + ip%HeIIcnt
  end if
#endif  

!---------------------------------------------------------
! optical depths

  ip%tauHI = ip%cdfac * ip%HIcnt * ip%sigmaHI
  ip%tausum = ip%tauHI
  if (ip%tauHI .LT. TAU_TOL) then
     ip%HItaufac = ip%tauHI
  else
     ip%HItaufac = 1.0 - exp(-ip%tauHI)
  end if
  ip%taufacsum = ip%HItaufac  

#ifdef incHe
  ip%tauHeI = ip%cdfac * ip%HeIcnt * ip%sigmaHeI
  ip%tauHeII = ip%cdfac * ip%HeIIcnt * ip%sigmaHeII
  ip%tausum = ip%tauHI + ip%tauHeI + ip%tauHeII 
  if (ip%tauHeI .LT. TAU_TOL) then
     ip%HeItaufac = ip%tauHeI
  else
     ip%HeItaufac = 1.0 - exp(-ip%tauHeI)
  end if
  if (ip%tauHeII .LT. TAU_TOL) then
     ip%HeIItaufac = ip%tauHeII
  else
     ip%HeIItaufac = 1.0 - exp(-ip%tauHeII)
  end if
  ip%taufacsum = ip%HItaufac + ip%HeItaufac + ip%HeIItaufac
#endif 

!---------------------------------------------------------
! deposit photons

  ip%fracabsorb = 1.0d0 - exp(-ip%tausum)

  ip%pdeps = ip%pflux * ip%fracabsorb


  ! if more photons going in then absorbers
  if (ip%pdeps >= ip%HIcnt) then
     ip%pdeps = ip%HIcnt
     ip%xHII = 1.0d0
     ip%xHI = 0.0d0
  else
     ip%xHII = ip%xHII + ip%pdeps / ip%Hcnt
     ip%xHI = 1.0d0 - ip%xHII
  end if

  scalls = scalls + 1
     
  
end subroutine recombeulerint




!> construct derivative matrices
!================================



!> this routine evolves the ionization fractions forward in time by dt
!======================================================================
subroutine take_euler_ion_step(xH,DH,dt,xHe,DHe)

  real(r8b), intent(inout) :: xH(2)
  real(r8b), intent(in) :: DH(2,2)
  real(r8b), intent(in) :: dt
  real(r8b), intent(inout), optional :: xHe(3)
  real(r8b), intent(in), optional :: DHe(3,3)

  real(r8b) :: MH(2,2)
  real(r8b) :: MHe(3,3)

  MH = invertDH(DH,dt)
  xH = matmul(MH,xH)

  if ( present(xHe) .and. present(DHe) ) then

     MHe = invertDHe(DHe,dt)
     xHe = matmul(MHe,xHe)

  end if
  
end subroutine take_euler_ion_step


!> inverts the Hydrogen derivative matrix (actually 1 + DH*dt)
!======================================================================
function invertDH(DH,dt) result(MH)

  real(r8b), parameter :: one = 1.0d0
  real(r8b) :: DH(2,2)
  real(r8b) :: dt
  real(r8b) :: MH(2,2)

  real(r8b) :: GG,RR,detH

  GG = DH(2,1)
  RR = DH(1,2)

  detH = one + ( RR + GG ) * dt

  MH(1,1) = one + RR * dt
  MH(2,1) = GG * dt
  MH(1,2) = RR * dt
  MH(2,2) = one + GG * dt

  MH = MH / detH

end function invertDH

!> inverts the Helium derivative matrix (actually 1 + DHe*dt)
!======================================================================
function invertDHe(DHe,dt) result(MHe)

  real(r8b), parameter :: one = 1.0d0
  real(r8b) :: DHe(3,3)
  real(r8b) :: dt
  real(r8b) :: MHe(3,3)

  real(r8b) :: GGI,GGII,RRII,RRIII,detHe

  GGI   = DHe(2,1)
  GGII  = DHe(3,2)
  RRII  = DHe(1,2)
  RRIII = DHe(2,3)

  detHe = one + (GGI + GGII + RRII + RRIII) * dt + &
          (GGI * GGII + GGI * RRIII + RRII * RRIII) * (dt*dt)

  MHe(1,1) = one + dt * (GGII + RRII + RRIII + dt * RRII * RRIII)  
  MHe(2,1) = dt * GGI * (one + dt * RRIII)
  MHe(3,1) = dt * dt * GGI * GGII
  
  MHe(1,2) = dt * RRII * (one + dt * RRIII)
  MHe(2,2) = (one + GGI * dt) * (one + dt * RRIII)
  MHe(3,2) = dt * (one + dt * GGI) * GGII

  MHe(1,3) = dt * dt * RRII * RRIII
  MHe(2,3) = dt * (one + dt * GGI) * RRIII
  MHe(3,3) = one + dt * (GGI + GGII + (dt * GGI * GGII) + RRII)

  MHe = MHe / detHe

end function invertDHe



end module euler_mod

