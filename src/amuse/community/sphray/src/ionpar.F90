!> \file ionpar.F90

!> \brief the module that handles ionization particles
!<

module ionpar_mod
use myf03_mod
use gadget_general_class
use physical_constants_mod, only: HI_th_Hz, HI_th_erg, HeI_th_erg, HeII_th_erg
use atomic_rates_mod, only: atomic_rates_type
use particle_system_mod, only: particle_type
use ray_mod
use raylist_mod
use b2cd_mod, only: b2cdfac
use cen_atomic_rates_mod, only: Verner_HI_photo_cs
use cen_atomic_rates_mod, only: Osterbrok_HeI_photo_cs
use cen_atomic_rates_mod, only: Osterbrok_HeII_photo_cs
use cen_atomic_rates_mod, only: Haiman_Bremss_cool
use cen_atomic_rates_mod, only: Haiman_Comp_Heol
use global_mod, only: GV
use hui_gnedin_atomic_rates_mod
implicit none


!> ionization particle type. 
!---------------------------
type ionpart_type

   ! quantities taken from standard particle

   real(r8b) :: pos(3)       !< x,y,z coordinates
   real(r8b) :: vel(3)       !< x,y,z velocities
   integer(i8b) :: id        !< particle id
   real(r8b) :: mass         !< particle mass
   real(r8b) :: T            !< temperature in K       
   real(r8b) :: rho          !< density
   real(r8b) :: hsml         !< smoothing length

   real(r8b) :: H_mf         !< Hydrogen mass fraction
   real(r8b) :: He_mf        !< Helium mass fraction

   real(r8b) :: xHI          !< nHI  / nH
   real(r8b) :: xHII         !< nHII / nH
   real(r8b) :: xHeI         !< nHeI   / nHe 
   real(r8b) :: xHeII        !< nHeII  / nHe 
   real(r8b) :: xHeIII       !< nHeIII / nHe

   integer(i8b) :: lasthit   !< last ray to cross this particle
   integer(i8b) :: index     !< index of the particle in the psys


   ! quantities initialized before solver is called 

   real(r8b) :: xHI_in       !< initial xHI
   real(r8b) :: xHII_in      !< initial xHII
   real(r8b) :: xHeI_in      !< initial xHeI
   real(r8b) :: xHeII_in     !< initial xHeII
   real(r8b) :: xHeIII_in    !< initial xHeIII
   real(r8b) :: T_in         !< initial T

   integer(i8b) :: iter      !< number of the iteration in the solver
   integer(i8b) :: rayn      !< ray number that impacted

   real(r8b) :: NeBckgnd     !< number density of metallic electrons
   real(r8b) :: Tcmb         !< background radiation field temperature

   real(r8b) :: gpercm3      !< density in cgs
   real(r8b) :: cm3          !< volume in cgs

   real(r8b) :: nH           !< number density of H
   real(r8b) :: nHe          !< number density of He
   real(r8b) :: Hcnt         !< number of H nuclei
   real(r8b) :: Hecnt        !< number of He nuclei

   real(r8b) :: dl           !< path length for particle
   logical   :: inside       !< is the source inside the particle? 

   ! ray/photo quantities initialied before solver is called

   integer(i8b) :: impact    !< index of the particle in the raylist
   real(r8b) :: d            !< distance along ray (code)
   real(r8b) :: b            !< impact parameter (code)
   real(r8b) :: bnorm        !< impact parameter normalized to smoothing length
   real(r8b) :: cdfac        !< column depth conversion factor

   real(r8b) :: dt_code      !< time step for this par (code)
   real(r8b) :: dt_s         !< time step for this par (s)

   real(r8b) :: penrg        !< energy of one photon in the ray (ergs)
   real(r8b) :: pflux        !< photons per second arriving at the particle
   real(r8b) :: pdeps        !< total photons deposited into particle
   real(r8b) :: pdeps_eq     !< photons deposited in equilibrium conditions

   real(r8b) :: DH(2,2)      !< dxH/dt array
   real(r8b) :: DHe(3,3)     !< dxHe/dt array


   real(r8b) :: sigmaHI      !< HI photo cross-section
   real(r8b) :: sigmaHeI     !< HeI photo cross-section
   real(r8b) :: sigmaHeII    !< HeII photo cross-section

   real(r8b) :: fracabsorb   !< fraction of all photons absorbed
   real(r8b) :: pdepr        !< photons per second deposited in particle


!   real(r8b) :: HIcolions   !< total HI collisional ionizations in the particle
!   real(r8b) :: HeIcolions  !< total HeI collisional ionizations in the particle
!   real(r8b) :: HeIIcolions !< total HeII collisional ionizations in the particle

!   real(r8b) :: HIIrcmbsB   !< total HII recombinations excluding to n=1 lvl
!   real(r8b) :: HIIrcmbsA   !< total HII recombinations to all levels
!   real(r8b) :: HeIIrcmbsB  !< total HeII recombinations excluding to n=1 lvl
!   real(r8b) :: HeIIrcmbsA  !< total HeII recombinations to all levels
!   real(r8b) :: HeIIIrcmbsB !< total HeIII recombinations excluding to n=1 lvl
!   real(r8b) :: HeIIIrcmbsA !< total HeIII recombinations to all levels



   real(r8b) :: ne           !< number density of electrons
   real(r8b) :: dnedt        !< time rate of change of ne

   real(r8b) :: nHI          !< number density of HI
   real(r8b) :: nHII         !< number density of HII
   real(r8b) :: nHeI         !< number density of HeI
   real(r8b) :: nHeII        !< number density of HeII
   real(r8b) :: nHeIII       !< number density of HeIII

   real(r8b) :: HIcnt        !< number of HI atoms
   real(r8b) :: HIIcnt       !< number of HII atoms
   real(r8b) :: HeIcnt       !< number of HeI atoms
   real(r8b) :: HeIIcnt      !< number of HeII atoms
   real(r8b) :: HeIIIcnt     !< number of HeIII atoms
   real(r8b) :: Allcnt       !< number of all photo absorbing species (HI,HeI,HeII)

   real(r8b) :: tauHI        !< HI optical depth
   real(r8b) :: tauHeI       !< HeI optical depth
   real(r8b) :: tauHeII      !< HeII optical depth
   real(r8b) :: tausum       !< sum of all optical depths

   real(r8b) :: HItaufac     !< 1-exp(-tauHI)
   real(r8b) :: HeItaufac    !< 1-exp(-tauHeI)
   real(r8b) :: HeIItaufac   !< 1-exp(-tauHeII)
   real(r8b) :: taufacsum    !< sum of all tau factors

   real(r8b) :: HIfrac       !< fraction of photons absorbed by HI
   real(r8b) :: HeIfrac      !< fraction of photons absorbed by HeI
   real(r8b) :: HeIIfrac     !< fraction of photons absorbed by HeII

   real(r8b) :: gammaHI      !< HI   photoionization rate
   real(r8b) :: gammaHeI     !< HeI  photoionization rate
   real(r8b) :: gammaHeII    !< HeII photoionization rate
   real(r8b) :: gammasum     !< sum of all photoionization rates

   real(r8b) :: GG           !< HIci * ne
   real(r8b) :: GGp          !< HI photoionization rate + HIci * ne
   real(r8b) :: RRa          !< HII recombination rate case A * ne
   real(r8b) :: RRb          !< HII recombination rate case B * ne

   real(r8b) :: GGI          !< HeIci * ne
   real(r8b) :: GGII         !< HeIIci * ne
   real(r8b) :: GGIp         !< HeI photoionization rate + HeIci * ne
   real(r8b) :: GGIIp        !< HeII photoionization rate + HeIIci * ne
   real(r8b) :: RRIIa        !< HeII recombination rate case A * ne
   real(r8b) :: RRIIIa       !< HeIII recombination rate case A * ne
   real(r8b) :: RRIIb        !< HeII recombination rate case B * ne
   real(r8b) :: RRIIIb       !< HeIII recombination rate case B * ne

   real(r8b) :: CIC          !< collisional ionization cooling rate
   real(r8b) :: CEC          !< collisional excitation cooling rate
   real(r8b) :: RCC          !< recombination cooling rate
   real(r8b) :: PH           !< photo heating rate
   real(r8b) :: BREM         !< bremsstrahlung cooling rate
   real(r8b) :: COMP         !< compton heating/cooling rate
   real(r8b) :: COOL         !< total cooling function w/o photo heating
   real(r8b) :: COOLp        !< total cooling function w photo heating

   real(r8b) :: tion         !< ionization time (s)
   real(r8b) :: tcool        !< cooling time (s)
   real(r8b) :: u            !< energy per unit mass of particle (ergs)
   real(r8b) :: dudt         !< time rate of change of energy
   real(r8b) :: dTdt         !< time rate of change of temperature
   
   character(200) :: strtag  !< labels code landmarks for debugging 

end type ionpart_type

 
!> particle type for solving the analytic ionization equations
!---------------------------------------------------------------
type bckgnd_particle_type
   real(r8b) :: H_mf
   real(r8b) :: T
   real(r8b) :: rho_cgs
   real(r8b) :: nH
   real(r8b) :: RC
   real(r8b) :: CI
   real(r8b) :: gammaHI
   real(r8b) :: y
   real(r8b) :: R
   real(r8b) :: Q
   real(r8b) :: P
   real(r8b) :: d
   real(r8b) :: xHI
   real(r8b) :: xHII
end type bckgnd_particle_type


!> particle type for calculating optical depth
!---------------------------------------------------------------
type tau_particle_type
   real(r8b) :: H_mf
   real(r8b) :: xHI
   real(r8b) :: xHII
   real(r8b) :: hsml
   real(r8b) :: mass_cgs
   real(r8b) :: Hcnt
   real(r8b) :: HIcnt
   real(r8b) :: d
   real(r8b) :: b
   real(r8b) :: bnorm
   real(r8b) :: cdfac
   real(r8b) :: tauHI_th
end type tau_particle_type




real(r8b), parameter, private :: zero = 0.0d0
real(r8b), parameter, private :: one = 1.0d0
real(r8b), parameter, private :: two = 2.0d0
real(r8b), parameter, private :: four = 4.0d0


contains


!> initializes a particle used just for calculating the optical depth
!-----------------------------------------------------------------------
function initialize_tau_particle( par, intersection ) result( tpar )
  type(particle_type), intent(in) :: par
  type(intersection_type), intent(in) :: intersection
  type(tau_particle_type) :: tpar
  
  type(gadget_constants_type) :: gconst
  logical, save :: first = .true.
  real(r8b), save :: sigmaHI_th

  if (first) then
     sigmaHI_th = Verner_HI_photo_cs(one)    
  endif

#ifdef incHmf
  tpar%H_mf   = par%Hmf
#else
  tpar%H_mf   = GV%H_mf
#endif

  tpar%xHI  = par%xHI
  tpar%xHII = par%xHII
  tpar%hsml = par%hsml

  tpar%mass_cgs = par%mass * GV%cgs_mass
  tpar%Hcnt     = tpar%mass_cgs * tpar%H_mf / gconst%PROTONMASS
  tpar%HIcnt    = tpar%Hcnt * tpar%xHI       

  tpar%d = intersection%d   ! distance along ray
  tpar%b = intersection%b   ! impact parameter
      
  tpar%bnorm = tpar%b / tpar%hsml
  tpar%cdfac = b2cdfac(tpar%bnorm, tpar%hsml, GV%cgs_len)

  if (tpar%xHI > zero) then
     tpar%tauHI_th = tpar%cdfac * tpar%HIcnt * sigmaHI_th
  else
     tpar%tauHI_th = zero
  end if
  
  first = .false.
  
end function initialize_tau_particle



!> initializes a particle used just for calculating the xH's
!-----------------------------------------------------------------------
function initialize_bckgnd_particle( par, gammaHI ) result( bpar )
  type(particle_type), intent(in) :: par
  real(r8b), intent(in) :: gammaHI
  type(bckgnd_particle_type) :: bpar
  type(gadget_constants_type) :: gconst

#ifdef incHmf
  bpar%H_mf   = par%Hmf
#else
  bpar%H_mf   = GV%H_mf
#endif

  bpar%gammaHI = gammaHI


  bpar%rho_cgs = par%rho * GV%cgs_rho
  bpar%nH      = bpar%rho_cgs * bpar%H_mf / gconst%PROTONMASS

  bpar%T = par%T
  bpar%RC   = Hui_HII_recombA( bpar%T )
  bpar%CI   = Hui_HI_col_ion( bpar%T )
  bpar%xHI  = par%xHI
  bpar%xHII = par%xHII
  
end function initialize_bckgnd_particle


!> implements analytic solution for xH's in background particles
!-----------------------------------------------------------------------
function set_bckgnd_particle_xH_eq( bpar ) result(err)
  type(bckgnd_particle_type) :: bpar
  integer(i4b) :: err

  bpar%y = 0.0d0

  bpar%R = -( bpar%CI + bpar%RC ) * bpar%nH
  bpar%Q = bpar%CI * bpar%nH - bpar%gammaHI - (bpar%CI + bpar%RC) * bpar%nH * bpar%y
  bpar%P = bpar%gammaHI + bpar%CI * bpar%nH * bpar%y

  bpar%d = bpar%Q * bpar%Q - four * bpar%R * bpar%P

  bpar%xHII = ( -bpar%Q - sqrt(bpar%D) ) / (two * bpar%R)
  bpar%xHI  = one - bpar%xHII

  if (bpar%xHII > one .or. bpar%xHII < zero) then
     err = -1
  else
     err = 0
  endif

end function set_bckgnd_particle_xH_eq



!> implements analytic solution for xH's in ion particles
!-----------------------------------------------------------------------
subroutine set_ionpar_xH_eq( ipar )
  type(ionpart_type), intent(inout) :: ipar  !< output ionization particle

  real(r8b) :: R, Q, P, D
  real(r8b) :: CI, RC
  real(r8b) :: y

  ipar%strtag = "in set_ionpar_xH_eq_enter"
  call check_x(ipar)

  y = 0.0d0
  CI = Hui_HI_col_ion( ipar%T )
  RC = Hui_HII_recombA( ipar%T )
  
  R = -( CI + RC ) * ipar%nH
  Q = CI * ipar%nH - ipar%gammaHI - (CI + RC) * ipar%nH * y
  P = ipar%gammaHI + CI * ipar%nH * y

  D = Q*Q - 4.0d0 * R * P

  ipar%xHII = ( - Q - sqrt(D)) / (2.0d0 * R)
  if (ipar%xHII > 1.0d0 .or. ipar%xHII < 0.0d0) then
     write(*,*) 'shit'
     stop
  endif

  ipar%xHI = 1.0d0 - ipar%xHII

  ipar%strtag = "in set_ionpar_xH_eq_exit"
  call check_x(ipar)

end subroutine set_ionpar_xH_eq




!> copies the basic particle data into an ionization particle
!-----------------------------------------------------------------------
subroutine par2ionpar(par,ipar,index)
 type(particle_type), intent(in) :: par     !< input particle
 type(ionpart_type), intent(inout) :: ipar  !< output ionization particle
 integer(i8b) :: index                      !< par = psys%par(index)

 ipar%pos     = par%pos

#ifdef incVel
 ipar%vel     = par%vel
#else
 ipar%vel     = 0.0
#endif

 ipar%id      = par%id
 ipar%mass    = par%mass
 ipar%T       = par%T
 ipar%rho     = par%rho

 ipar%xHI     = par%xHI
 ipar%xHII    = par%xHII

 ipar%hsml    = par%hsml

#ifdef incHmf
  ipar%H_mf   = par%Hmf
#else
  ipar%H_mf   = GV%H_mf
#endif

#ifdef incHemf
  ipar%He_mf  = par%Hemf
#else
  ipar%He_mf  = GV%He_mf
#endif


#ifdef incHe
 ipar%xHeI    = par%xHeI
 ipar%xHeII   = par%xHeII
 ipar%xHeIII  = par%xHeIII
#else
 ipar%xHeI    = 0.0d0
 ipar%xHeII   = 0.0d0
 ipar%xHeIII  = 0.0d0
#endif

 ipar%lasthit = par%lasthit
 ipar%index   = index

end subroutine par2ionpar


!> copies the ionization particle data into a basic particle 
!-----------------------------------------------------------------------
subroutine ionpar2par(ipar,par)
 type(ionpart_type), intent(in) :: ipar    !< input ionization particle
 type(particle_type), intent(inout) :: par !< output particle

 par%T       = ipar%T

 par%xHI     = ipar%xHI
 par%xHII    = ipar%xHII

#ifdef incHe
 par%xHeI    = ipar%xHeI
 par%xHeII   = ipar%xHeII
 par%xHeIII  = ipar%xHeIII
#endif

 par%lasthit = ipar%lasthit

end subroutine ionpar2par


!> initializes the ionization particle values
!-----------------------------------------------------------------------
subroutine initialize_ionpar(ipar,par,index,srcray,He,raylist,impact)

  type(ionpart_type), intent(inout) :: ipar           !< ionization particle
  type(particle_type), intent(in) :: par              !< standard particle
  integer(i8b) :: index                               !< par = psys%par(index)
  logical, intent(in) :: srcray                       !< source ray update ?
  logical, intent(in) :: He                           !< update Helium ?
  type(raylist_type), intent(in), optional :: raylist !< optional raylist
  integer(i8b), intent(in), optional :: impact        !< optional impact number

  type(gadget_constants_type) :: gconst  
  real(r8b) :: mass_cgs
  real(r8b) :: rho_cgs


  call par2ionpar(par,ipar,index)

  ! store initial values
  !-------------------------
  ipar%xHI_in  = ipar%xHI
  ipar%xHII_in = ipar%xHII
  if (He) then
     ipar%xHeI_in   = ipar%xHeI
     ipar%xHeII_in  = ipar%xHeII
     ipar%xHeIII_in = ipar%xHeIII
  else
     ipar%xHeI_in   = 0.0d0
     ipar%xHeII_in  = 0.0d0
     ipar%xHeIII_in = 0.0d0
  end if
  ipar%T_in = ipar%T

  ! set values that are static during the update
  !-----------------------------------------------
  ipar%rayn     = GV%rayn

  ipar%NeBckgnd = GV%NeBackground
  ipar%Tcmb     = GV%Tcmb_cur

  mass_cgs      = ipar%mass * GV%cgs_mass 
  rho_cgs       = ipar%rho *  GV%cgs_rho 

  ipar%gpercm3  = rho_cgs
  ipar%cm3      = mass_cgs / rho_cgs

  ipar%nH       = rho_cgs  * ipar%H_mf / gconst%PROTONMASS
  ipar%Hcnt     = mass_cgs * ipar%H_mf / gconst%PROTONMASS

  if (He) then
     ipar%nHe   = rho_cgs  * ipar%He_mf / (4 * gconst%PROTONMASS)
     ipar%Hecnt = mass_cgs * ipar%He_mf / (4 * gconst%PROTONMASS)
  else
     ipar%nHe   = 0.0d0
     ipar%Hecnt = 0.0d0
  end if

  ! initialize solver iteration counter to zero
  !---------------------------------------------
  ipar%iter = 0


  ! set ray specific quantities 
  !---------------------------------------------
  if (present(raylist)) then

     if (.not. present(impact)) stop "raylist but no impact in initialize_ionpar"
     ipar%impact = impact  ! which impact in the raylist


     ipar%d = raylist%intersection(impact)%d   ! distance along ray
     ipar%b = raylist%intersection(impact)%b   ! impact parameter
 
!!$     if ( sqrt( sum( (raylist%ray%start - ipar%pos)**2 ) ) <= ipar%hsml ) then
!!$        ipar%inside=.true.
!!$     else
!!$        ipar%inside=.false.
!!$     endif
!!$
!!$     if (raylist%nnb == 1) then 
!!$        ipar%dl = 2 * ipar%hsml
!!$
!!$     else 
!!$        if (impact == 1) then 
!!$           if (ipar%inside) then
!!$              ipar%dl = ipar%hsml
!!$           else
!!$              ipar%dl = 0.5 * (raylist%intersection(2)%d + ipar%d )  
!!$!              ipar%dl = raylist%intersection(2)%d - ipar%d  
!!$!              ipar%dl = 2 * ipar%hsml
!!$           endif
!!$        else if (impact == raylist%nnb) then
!!$           ipar%dl = ipar%d - raylist%intersection(raylist%nnb-1)%d
!!$        else
!!$           ipar%dl = 0.5 * ( raylist%intersection(impact+1)%d - &
!!$                             raylist%intersection(impact-1)%d )
!!$        endif
!!$     endif

     
     ipar%bnorm = ipar%b / ipar%hsml
     ipar%cdfac = b2cdfac(ipar%bnorm,ipar%hsml,GV%cgs_len)

     ipar%dt_code = (GV%itime - ipar%lasthit) * GV%dt_code
     ipar%dt_s    = (GV%itime - ipar%lasthit) * GV%dt_s
     
     ipar%pdeps    = 0.0d0
     ipar%pdeps_eq = 0.0d0

     ipar%penrg = raylist%ray%enrg

     ipar%sigmaHI = Verner_HI_photo_cs(raylist%ray%freq)    
     if (He) then
        ipar%sigmaHeI = Osterbrok_HeI_photo_cs(raylist%ray%freq * HI_th_Hz)    
        ipar%sigmaHeII = Osterbrok_HeII_photo_cs(raylist%ray%freq * HI_th_Hz)
     else
        ipar%sigmaHeI = 0.0d0
        ipar%sigmaHeII = 0.0d0
     end if
     
     if (srcray) then        
        ipar%pflux = raylist%ray%pcnt / ipar%dt_s 
     else
        ipar%pflux = raylist%ray%pcnt 
     end if

  end if


  ipar%strtag = "in initialize_ionpar"
  call check_x(ipar)
     
end subroutine initialize_ionpar






!> initializes a non photo update particle
!-----------------------------------------------------------------------
subroutine initialize_non_photo_ionpar(ipar,par,index,He)

  type(ionpart_type), intent(inout) :: ipar           !< ionization particle
  type(particle_type), intent(in) :: par              !< standard particle
  integer(i8b) :: index                               !< par = psys%par(index)
  logical, intent(in) :: He                           !< update Helium ?

  type(gadget_constants_type) :: gconst  
  real(r8b) :: mass_cgs
  real(r8b) :: rho_cgs


  call par2ionpar(par,ipar,index)

  ! store initial values
  !-------------------------
  ipar%xHI_in  = ipar%xHI
  ipar%xHII_in = ipar%xHII
  if (He) then
     ipar%xHeI_in   = ipar%xHeI
     ipar%xHeII_in  = ipar%xHeII
     ipar%xHeIII_in = ipar%xHeIII
  else
     ipar%xHeI_in   = 0.0d0
     ipar%xHeII_in  = 0.0d0
     ipar%xHeIII_in = 0.0d0
  end if
  ipar%T_in = ipar%T

  ! set values that are static during the update
  !-----------------------------------------------

  ipar%NeBckgnd = GV%NeBackground
  ipar%Tcmb     = GV%Tcmb_cur

  mass_cgs      = ipar%mass * GV%cgs_mass 
  rho_cgs       = ipar%rho *  GV%cgs_rho 

  ipar%gpercm3  = rho_cgs
  ipar%cm3      = mass_cgs / rho_cgs

  ipar%nH       = rho_cgs  * ipar%H_mf / gconst%PROTONMASS
  ipar%Hcnt     = mass_cgs * ipar%H_mf / gconst%PROTONMASS

  if (He) then
     ipar%nHe   = rho_cgs  * ipar%He_mf / (4 * gconst%PROTONMASS)
     ipar%Hecnt = mass_cgs * ipar%He_mf / (4 * gconst%PROTONMASS)
  else
     ipar%nHe   = 0.0d0
     ipar%Hecnt = 0.0d0
  end if

  ! initialize solver iteration counter to zero
  !---------------------------------------------
  ipar%iter = 0

  ipar%dt_code = (GV%itime - ipar%lasthit) * GV%dt_code
  ipar%dt_s    = (GV%itime - ipar%lasthit) * GV%dt_s
     
  ipar%strtag = "in initialize_non_photo_ionpar"
  call check_x(ipar)
     
end subroutine initialize_non_photo_ionpar









!> sets the photoionization rate for an ionization particle
!================================================================
subroutine set_taus(ip,He)

  real(r8b), parameter :: zero = 0.0d0
  real(r8b), parameter :: one = 1.0d0
  real(r8b), parameter :: TAU_LOW = 1.0d-4
  real(r8b), parameter :: TAU_HIGH = 3.0d1
  type(ionpart_type), intent(inout) :: ip !< ionization particle  
  logical, intent(in) :: He

  logical :: HI,HeI,HeII


  ! check which species will be absorbing
  !---------------------------------------
  HI   = .false.
  HeI  = .false.
  HeII = .false.

  if (ip%xHI > zero) HI = .true.

  if (He) then
     if (ip%xHeI  > zero .and. ip%penrg > HeI_th_erg)  HeI = .true.
     if (ip%xHeII > zero .and. ip%penrg > HeII_th_erg) HeII = .true.
  end if


  ! calculate atom counts
  !---------------------------------------
  ip%HIcnt = ip%Hcnt * ip%xHI       
  if (He) then
     ip%HeIcnt  = ip%Hecnt * ip%xHeI        
     ip%HeIIcnt = ip%Hecnt * ip%xHeII          
  else
     ip%HeIcnt = zero
     ip%HeIIcnt = zero
  end if

  ip%Allcnt = ip%HIcnt + ip%HeIcnt + ip%HeIIcnt


  ! calculate taus
  !---------------------------------------
  if (HI) then
     ip%tauHI = ip%cdfac * ip%HIcnt * ip%sigmaHI
!     write(*,*) "tau int: ", ip%tauHI
!     ip%tauHI = ip%dl * GV%cgs_len * ip%nH * ip%xHI * ip%sigmaHI
!     write(*,*) "tau path:", ip%tauHI
!     write(*,*) 
  else
     ip%tauHI = zero
  end if

  !---------------------------------------
  if (HeI) then
     ip%tauHeI = ip%cdfac * ip%HeIcnt * ip%sigmaHeI
  else
     ip%tauHeI = zero
  end if
     
  !---------------------------------------
  if (HeII) then
     ip%tauHeII = ip%cdfac * ip%HeIIcnt * ip%sigmaHeII
  else 
     ip%tauHeII = zero
  end if
  
  ip%tausum = ip%tauHI + ip%tauHeI + ip%tauHeII


  ! calculate absorption ratios
  !--------------------------------
  if (HI) then

     if (ip%tauHI < TAU_LOW) then
        ip%HItaufac = ip%tauHI
     else if (ip%tauHI > TAU_HIGH) then
        ip%HItaufac = one
     else
        ip%HItaufac = one - exp(-ip%tauHI)
     end if

  else
     
     ip%HItaufac = zero

  end if

  !-------------------------------------
  if (HeI) then

     if (ip%tauHeI < TAU_LOW) then
        ip%HeItaufac = ip%tauHeI
     else if (ip%tauHeI > TAU_HIGH) then
        ip%HeItaufac = one
     else
        ip%HeItaufac = one - exp(-ip%tauHeI)
     end if

  else
     
     ip%HeItaufac = zero

  end if

  !-------------------------------------
  if (HeII) then

     if (ip%tauHeII < TAU_LOW) then
        ip%HeIItaufac = ip%tauHeII
     else if (ip%tauHeII > TAU_HIGH) then
        ip%HeIItaufac = one
     else
        ip%HeIItaufac = one - exp(-ip%tauHeII)
     end if

  else
     
     ip%HeIItaufac = zero

  end if

  ip%taufacsum = ip%HItaufac + ip%HeItaufac + ip%HeIItaufac

  if (ip%taufacsum > zero) then
  
     if (.not. He) then
        ip%HIfrac   = one
        ip%HeIfrac  = zero
        ip%HeIIfrac = zero
     else
        ip%HIfrac   = ip%HItaufac   / ip%taufacsum
        ip%HeIfrac  = ip%HeItaufac  / ip%taufacsum
        ip%HeIIfrac = ip%HeIItaufac / ip%taufacsum
     end if

  else

     ip%HIfrac   = zero
     ip%HeIfrac  = zero
     ip%HeIIfrac = zero
     
  end if

end subroutine set_taus


!> sets the photoionization rate for an ionization particle
!================================================================
subroutine set_gammas(ip,He)

  real(r8b), parameter :: zero = 0.0d0
  real(r8b), parameter :: one = 1.0d0
  real(r8b), parameter :: TAU_LOW = 1.0d-4
  real(r8b), parameter :: TAU_HIGH = 3.0d+1
  type(ionpart_type), intent(inout) :: ip !< ionization particle  
  logical, intent(in) :: He

  
  call set_taus(ip,He)

  ! photoionization rates
  !------------------------
  if (ip%tausum > zero) then

     if (ip%tausum < TAU_LOW) then
        ip%fracabsorb = ip%tausum
     else if (ip%tausum > TAU_HIGH) then
        ip%fracabsorb = one
     else
        ip%fracabsorb = one - exp(-ip%tausum)
     end if

     ip%pdepr    = ip%pflux * ip%fracabsorb
     ip%gammasum = ip%pdepr / ip%Allcnt

     ip%gammaHI = ip%gammasum * ip%HIfrac 
     if (He) then
        ip%gammaHeI  = ip%gammasum * ip%HeIfrac 
        ip%gammaHeII = ip%gammasum * ip%HeIIfrac 
     end if

  else

     ip%gammasum = zero
     ip%gammaHI = zero
     if (He) then
        ip%gammaHeI = zero
        ip%gammaHeII = zero
     end if
     return

  end if
  
end subroutine set_gammas




!> sets the time rate of change of the electron number density
!================================================================
subroutine set_dnedt(ip,photo,caseA,He)

  type(ionpart_type), intent(inout) :: ip   !< ionization particle
  logical, intent(in) :: photo
  logical, intent(in) :: caseA(2)
  logical, intent(in) :: He

  real(r8b) :: GG,RR
  real(r8b) :: GGI,GGII,RRII,RRIII
  
  if (photo) then
     GG = ip%GGp
  else
     GG = ip%GG
  end if

  if (caseA(1)) then
     RR = ip%RRa
  else
     RR = ip%RRb
  end if

  ip%dnedt = GG * ip%nHI - RR * ip%nHII

  if (He) then

     if (photo) then
        GGI  = ip%GGIp
        GGII = ip%GGIIp
     else
        GGI  = ip%GGI
        GGII = ip%GGII
     end if

     if (caseA(2)) then
        RRII  = ip%RRIIa
        RRIII = ip%RRIIIa
     else
        RRII  = ip%RRIIb
        RRIII = ip%RRIIIb
     end if

     ip%dnedt = ip%dnedt +  (GGI   * ip%nHeI ) + (GGII  * ip%nHeII )   &
                         -  (RRII  * ip%nHeII) - (RRIII * ip%nHeIII)

  end if

end subroutine set_dnedt
 
!> sets the GG's and RR's from the atomic rates, ne's and gammas
!=================================================================
subroutine set_ionization_func(ip,k,photo,caseA,He)

  type(ionpart_type), intent(inout) :: ip   !< ionization particle  
  type(atomic_rates_type), intent(in) :: k  !< rates
  logical, intent(in) :: photo
  logical, intent(in) :: caseA(2)
  logical, intent(in) :: He
 
  ip%GG  = k%HIci * ip%ne 
  if (photo) ip%GGp = ip%GG + ip%gammaHI
  
  ip%RRb = k%HIIrcB * ip%ne
  if (caseA(1)) ip%RRa = k%HIIrcA * ip%ne

  if (He) then
     ip%GGI    = k%HeIci  * ip%ne 
     ip%GGII   = k%HeIIci * ip%ne 
     if (photo) then
        ip%GGIp  = ip%GGI  + ip%gammaHeI  
        ip%GGIIp = ip%GGII + ip%gammaHeII  
     end if

     ip%RRIIb  = k%HeIIrcB  * ip%ne
     ip%RRIIIb = k%HeIIIrcB * ip%ne

     if (caseA(2)) then
        ip%RRIIa  = k%HeIIrcA  * ip%ne
        ip%RRIIIa = k%HeIIIrcA * ip%ne
     end if
  end if

end subroutine set_ionization_func

 
!> sets the global cooling rate for an ionization particle
!================================================================
subroutine set_cooling_func(ip,k,photo,caseA,He)

  type(ionpart_type), intent(inout) :: ip   !< ionization particle  
  type(atomic_rates_type), intent(in) :: k  !< rates
  logical, intent(in) :: photo
  logical, intent(in) :: caseA(2)
  logical, intent(in) :: He

  ! CIC  = collisional ionization cooling
  ! CEC  = collisional excitation cooling
  ! RCC  = recombination cooling 
  ! BREM = free free cooling
  ! COMP = compton cooling
  ! PH   = photo heating  

  ! cooling from Hydrogen
  !-----------------------
  ip%PH  = 0.0d0
  ip%CIC = ( k%HIcic * ip%nHI  ) * ip%ne 
  ip%CEC = ( k%HIcec * ip%nHI  ) * ip%ne

  if (caseA(1)) then
     ip%RCC = ( k%HIIrccA * ip%nHII ) * ip%ne
  else
     ip%RCC = ( k%HIIrccB * ip%nHII ) * ip%ne
  end if

  if (photo) then
     ip%PH = ip%pdepr * ip%HIfrac * (ip%penrg - HI_th_erg)
  end if

  ! cooling from Helium (yes there is an extra ne in He CEC)
  !----------------------------------------------------------
  if (He) then
     ip%CIC = ip%CIC + ( k%HeIcic  * ip%nHeI  + k%HeIIcic * ip%nHeII ) * ip%ne
     ip%CEC = ip%CEC + ( k%HeIcec  * ip%nHeII * ip%ne + k%HeIIcec * ip%nHeII ) * ip%ne

     if (caseA(2)) then
        ip%RCC = ip%RCC + ( k%HeIIrccA * ip%nHeII + k%HeIIIrccA * ip%nHeIII ) * ip%ne 
     else
        ip%RCC = ip%RCC + ( k%HeIIrccB * ip%nHeII + k%HeIIIrccB * ip%nHeIII ) * ip%ne 
     end if

     if (photo) then
        ip%PH = ip%PH + ip%pdepr * ip%HeIfrac  * (ip%penrg - HeI_th_erg ) 
        ip%PH = ip%PH + ip%pdepr * ip%HeIIfrac * (ip%penrg - HeII_th_erg) 
     end if
  end if

  ip%PH = ip%PH / ip%cm3
  ip%BREM = Haiman_Bremss_cool(ip%T,ip%nHII,ip%nHeII,ip%nHeIII,ip%ne)
  ip%COMP = Haiman_Comp_Heol(ip%T,ip%Tcmb,ip%ne)

  ip%COOL  = -( ip%CIC + ip%CEC + ip%RCC + ip%BREM + ip%COMP )
  ip%COOLp = ip%PH + ip%COOL 
 
  
end subroutine set_cooling_func 

!> sets the Hydrogen derivative matrix from the values in GGRR
!======================================================================
subroutine setDH(ip,photo,caseA)

  type(ionpart_type), intent(inout) :: ip   !< ionization particle  
  logical, intent(in) :: photo
  logical, intent(in) :: caseA
  real(r8b) :: GG,RR
 
  if (photo) then
     GG = ip%GGp
  else
     GG = ip%GG
  end if

  if (caseA) then
     RR = ip%RRa
  else
     RR = ip%RRb
  end if

  ip%DH(1,1) = -GG
  ip%DH(2,1) = GG

  ip%DH(1,2) = RR
  ip%DH(2,2) = -RR

end subroutine setDH


!> sets the Helium derivative matrix from the values in GGRR
!======================================================================
subroutine setDHe(ip,photo,caseA)

  real(r8b), parameter :: zero = 0.0d0
  type(ionpart_type), intent(inout) :: ip   !< ionization particle  
  logical, intent(in) :: photo
  logical, intent(in) :: caseA

  real(r8b) :: GGI,GGII,RRII,RRIII

  if (photo) then
     GGI  = ip%GGIp
     GGII = ip%GGIIp
  else
     GGI  = ip%GGI
     GGII = ip%GGII     
  end if

  if (caseA) then
     RRII  = ip%RRIIa
     RRIII = ip%RRIIIa
  else
     RRII  = ip%RRIIb
     RRIII = ip%RRIIIb
  end if

  ip%DHe(1,1) = -GGI
  ip%DHe(2,1) = GGI
  ip%DHe(3,1) = zero

  ip%DHe(1,2) = RRII
  ip%DHe(2,2) = -(GGII + RRII)
  ip%DHe(3,2) = GGII

  ip%DHe(1,3) = zero
  ip%DHe(2,3) = RRIII
  ip%DHe(3,3) = -RRII

  
end subroutine setDHe



!> prints ionization particle information to screen
!-----------------------------------------------------
subroutine ionpar2screen(ipar)

  type(ionpart_type), intent(inout) :: ipar  !< ionization particle
  
   95 format(A,T25,I15)
   100 format(A,T25,3F12.4)
   102 format(A,T25,2F12.4)
   105 format(A,T25,1F12.4)
   110 format(A,T25,3ES18.9)
   115 format(A,T25,2ES18.9)
   write(*,*) 
   write(*,*) "ID", ipar%id
   write(*,*) "pos", ipar%pos
   write(*,*) "vel", ipar%vel
   write(*,*) "hsml", ipar%hsml
   write(*,*) "rho", ipar%rho
   write(*,*) "mass", ipar%mass
   write(*,*) "T", ipar%T
   write(*,*) "lasthit", ipar%lasthit
   write(*,*) 
   write(*,*) "ray num", ipar%rayn
   write(*,*) "psys index", ipar%index
   write(*,*) "impact", ipar%impact
   write(*,*) "iteration", ipar%iter
   write(*,*) "raydist", ipar%d
   write(*,*) 
   write(*,*) "string tag = ", trim(ipar%strtag)
   write(*,*) "xHI+xHII", ipar%xHI + ipar%xHII
   write(*,*) "xHeI+xHeII+xHeIII", ipar%xHeI + ipar%xHeII + ipar%xHeIII
   write(*,'(A,T25,ES18.9)') "T_in", ipar%T_in
   write(*,110) "b, hsml, bnorm",ipar%b, ipar%hsml, ipar%bnorm
   write(*,*) 
   write(*,*) "penrg/(HIth,HeIth,HeIIth)"
   write(*,*) ipar%penrg / HI_th_erg, ipar%penrg / HeI_th_erg, &
              ipar%penrg / HeII_th_erg
   write(*,*) 
   write(*,*) "input(xHI,xHII / xHeI,xHeII,xHeIII)"
   write(*,*) ipar%xHI_in, ipar%xHII_in
   write(*,*) ipar%xHeI_in, ipar%xHeII_in, ipar%xHeIII_in
   write(*,*) 
   write(*,*) "current(xHI,xHII / xHeI,xHeII,xHeIII)"
   write(*,*) ipar%xHI, ipar%xHII
   write(*,*) ipar%xHeI, ipar%xHeII, ipar%xHeIII
   write(*,*) 
   write(*,*) "current(nHI,nHII / nHeI,nHeII,nHeIII)"
   write(*,*) ipar%nHI, ipar%nHII
   write(*,*) ipar%nHeI, ipar%nHeII, ipar%nHeIII
   write(*,*) 
   write(*,*) "optical depths / tau facs (HI,HeI,HeII)"
   write(*,*) ipar%tauHI, ipar%tauHeI, ipar%tauHeII
   write(*,*) ipar%HItaufac, ipar%HeItaufac, ipar%HeIItaufac
   write(*,*) 
   write(*,*) "nH,nHe,ne / Hcnt,Hecnt,dnedt"
   write(*,*)  ipar%nH, ipar%nHe, ipar%ne
   write(*,*)  ipar%Hcnt, ipar%Hecnt, ipar%dnedt
   write(*,*) 
   write(*,*) "HIcnt,HeIcnt,HeIIcnt / HIfrac,HeIfrac,HeIIfrac"
   write(*,*) ipar%HIcnt, ipar%HeIcnt, ipar%HeIIcnt 
   write(*,*) ipar%HIfrac,ipar%HeIfrac,ipar%HeIIfrac
   write(*,*)
   write(*,*) "sigmas (HI,HeI,HeII) / Gammas (HI,HeI,HeII)"
   write(*,*) ipar%sigmaHI, ipar%sigmaHeI, ipar%sigmaHeII
   write(*,*) ipar%gammaHI, ipar%gammaHeI, ipar%gammaHeII
   write(*,*)
   write(*,*) "GGp (HI,HeI,HeII) / RRb (HII,HeII,HeIII)"
   write(*,*) ipar%GGp, ipar%GGIp, ipar%GGIIp
   write(*,*) ipar%RRb, ipar%RRIIb, ipar%RRIIIb
   write(*,*)
   write(*,*) "CIC,CEC,RCC / PH,BREM,COMP / COOLp, COOL"
   write(*,*) ipar%CIC, ipar%CEC, ipar%RCC
   write(*,*) ipar%PH, ipar%BREM, ipar%COMP
   write(*,*) ipar%COOLp, ipar%COOL
   write(*,*) 
   write(*,*) "pflux, fracabsorb, pdepr"
   write(*,*) ipar%pflux, ipar%fracabsorb, ipar%pdepr
   write(*,*) 
   write(*,110) "pdepr,penrg,pdeps", ipar%pdepr, ipar%penrg, ipar%pdeps

   write(*,*)
   write(*,115) "cdfac,d", ipar%cdfac,  ipar%d
   write(*,*) 
   write(*,115) "gpercm3,cm3", ipar%gpercm3, ipar%cm3

   write(*,*) 
   write(*,115) "u,dudt", ipar%u, ipar%dudt 
   write(*,115) "t,dTdt", ipar%T, ipar%dTdt
   write(*,*) 
   write(*,115) "Allcnt,GammaSum", ipar%Allcnt,ipar%gammasum
   write(*,115) "tausum,taufacsum",ipar%tausum, ipar%taufacsum
   write(*,*) 
   write(*,115) "tion,tcool", ipar%tion, ipar%tcool
   write(*,115) "dt_code, dt_s", ipar%dt_code, ipar%dt_s

   write(*,'(A,T25,ES18.9)') "pdeps_eq", ipar%pdeps_eq

   ipar%strtag = ""

 end subroutine ionpar2screen


!================================================================
!> dummy checks the bounds on HII, HeII, HeIII, and T
subroutine check_x(ip) 

  real(r8b), parameter :: zero = 0.0d0
  real(r8b), parameter :: one = 1.0d0

  real(r8b), parameter :: TOL = 0.0d+0
  type(ionpart_type), intent(inout) :: ip  !< ionization particle
  logical :: bad

  100 format(A,I2,A)
  101 format(A,I2,A,F15.8)
  102 format(A,I15)
  103 format(A,2F15.8)

  bad = .false.
 
  if ( ip%xHI .LT. zero - TOL ) then
     bad = .true.
     write(*,100) "xHI < zero in check_x"
  end if
     
  if ( ip%xHI > one + TOL ) then
     bad = .true.
     write(*,*) "xHI = ", ip%xHI
     write(*,100) "xHI > one in check_x"
  end if

  if ( ip%xHII .LT. zero - TOL ) then
     bad = .true.
     write(*,100) "xHII < zero in check_x"
  end if
     
  if ( ip%xHII .GT. one + TOL ) then
     bad = .true.
     write(*,100) "xHII > one in check_x"
  end if
     

  if ( ip%xHeI .LT. zero - TOL ) then
     bad = .true.
     write(*,100) "xHeI < zero in check_x"
  end if
     
  if ( ip%xHeI .GT. one + TOL ) then
     bad = .true.
     write(*,100) "xHeI > one in check_x"
  end if

  if ( ip%xHeII .LT. zero - TOL ) then
     bad = .true.
     write(*,100) "xHeII < zero in check_x"
  end if
     
  if ( ip%xHeII .GT. one + TOL ) then
     bad = .true.
     write(*,100) "xHeII > one in check_x"
  end if

  if ( ip%xHeIII .LT. zero - TOL ) then
     bad = .true.
     write(*,100) "xHeIII < zero in check_x"
  end if
     
  if ( ip%xHeIII .GT. one + TOL ) then
     bad = .true.
     write(*,100) "xHeIII > one in check_x"
  end if

  if ( ip%T .LE. zero ) then
     bad = .true.
     write(*,*) "T <= zero in check_x"
  end if


  if (bad) then
     call ionpar2screen(ip)
     stop
  end if

end subroutine check_x


end module ionpar_mod
