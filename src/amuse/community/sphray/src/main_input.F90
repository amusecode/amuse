!> \file main_input.F90

!> \brief The module that calls the specific input routines
!<

module main_input_mod
use myf03_mod
use gadget_general_class
use gadget_public_input_mod
use gadget_cosmoBH_input_mod
use gadget_owls_input_mod
use gadget_vbromm_input_mod
use gadget_public_input_hdf5_mod

use update_particles_mod
use source_input_mod
use particle_system_mod

use atomic_rates_mod, only: get_atomic_rates
use global_mod, only: PLAN, GV, rtable, cmbT_k
use global_mod, only: psys, saved_gheads
use global_mod, only: set_dt_from_dtcode
implicit none

contains


!> Read in planning data from the header of all snapshots 
!========================================================
subroutine get_planning_data()  
  character(clen), parameter :: myname="get_planning_data"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=1

  call mywrite("getting planning data:", verb)
  call mywrite("",verb) 


  GV%Nsnaps = GV%EndSnapNum - GV%StartSnapNum + 1
  if (allocated(PLAN%snap)) deallocate(PLAN%snap)
  allocate( PLAN%snap(GV%StartSnapNum : GV%EndSnapNum) )

  ! branch on input type 
  !-------------------------------------------------------------------
  select case (GV%InputType)
  case(1)
     call get_planning_data_gadget_public()
  case(2)
     call get_planning_data_gadget_cosmoBH()
  case(3)
     call get_planning_data_gadget_owls()
  case(4)
     call get_planning_data_gadget_vbromm()
  case(5)
     call get_planning_data_gadget_public_hdf5()
  end select

  call get_planning_data_sources()

end subroutine get_planning_data


!> read in particle, box, and source data 
!============================================
subroutine readin_snapshot()
  character(clen), parameter :: myname="readin_snapshot"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt

  
  logical :: first
  real(r8b) :: MB
  integer(i8b) :: i
  real(r8b) :: a     !< scale factor
  real(r8b) :: h     !< Hubble paraemter (little H)
  real(r8b) :: Flux  !< photons / s from planar sources
  character(clen) :: snpbase
  character(clen) :: srcbase
  type(gadget_constants_type) :: gconst
  
  call mywrite("reading in particle and source snapshots:", verb-1)
  call mywrite("",verb-1) 

  ! set local variables
  !======================
  a = PLAN%snap(GV%CurSnapNum)%ScalefacAt 
  h = GV%LittleH

  ! report readin type
  !======================
  call mywrite('   input type = ', verb, adv=.false.)
  if (GV%InputType==1) then
     call mywrite(" Gadget-2 Public (SnapFormat=1)", verb)
  else if (GV%InputType==2) then
     call mywrite(" Gadget CosmoBH w/ ye and xHI (SnapFormat=1)", verb)
  else if (GV%InputType==3) then
     call mywrite(" Gadget OWLS/GIMIC HDF5", verb)
  else if (GV%InputType==4) then
     call mywrite(" Gadget Volker Bromm", verb)
  else if (GV%InputType==5) then
     call mywrite(" Gadget-2 Public HDF5", verb)
  end if
 
  ! read in the particle data
  !=============================
  first = .false.
  if (GV%CurSnapNum == GV%StartSnapNum) then
     first = .true.
  else
     GV%CurSnapNum = GV%CurSnapNum + 1
  end if

  ! gadget public (no ionization fractions)
  !---------------------------------------------------------------
  if (GV%InputType == 1) then

     if (first) then
        call read_Gpublic_particles()
        psys%par(:)%lasthit = 0
     else
        call update_particles()
     end if

  ! gadget cosmoBH w/ cooling (i.e. ye and xHI)
  !---------------------------------------------------------------
  else if (GV%InputType == 2) then

     if (first) then
        call read_GcosmoBH_particles()
        psys%par(:)%lasthit = 0
     else
        call update_particles()
     end if
     
  ! gadget OWLS/GIMIC HDF5
  !---------------------------------------------------------------
  else if (GV%InputType == 3) then 
     
     if (first) then
        call read_Gowls_particles()
        psys%par(:)%lasthit = 0
     else
        call update_particles()
     end if
     
  ! gadget w/ ions from Volker Bromm's group 
  !---------------------------------------------------------------
  else if (GV%InputType == 4) then

     if (first) then
        call read_Gvbromm_particles()
        psys%par(:)%lasthit = 0
     else
        call update_particles()
     end if

  ! gadget public HDF5 (no ionization fractions)
  !---------------------------------------------------------------
  else if (GV%InputType == 5) then

     if (first) then
        call read_Gpubhdf5_particles()
        psys%par(:)%lasthit = 0
     else
        call update_particles()
     end if

  ! not recognized
  !---------------------------------------------------------------
  else
     write(str,*) "input type, ", GV%InputType, "not recognized" 
     call myerr(str,myname,crash)
  end if


  ! copy over box properties
  !==================================================
  psys%box%tops = GV%BoxUprs
  psys%box%bots = GV%BoxLwrs

  psys%box%lens    = GV%BoxUprs - GV%BoxLwrs
  psys%box%lens_cm = psys%box%lens * GV%cgs_len

  psys%box%vol    = product( psys%box%lens )
  psys%box%vol_cm = product( psys%box%lens_cm )

  psys%box%tbound = GV%BndryCond
  psys%box%bbound = GV%BndryCond


  ! read in the source data
  !============================================
  call read_src_snapshot()
  call order_sources_lum(psys%src)
  psys%src%lastemit = GV%rayn

  
  ! these quantities track the photoionization rate.  they are 
  ! rezeroed at inputs (because new source files are loaded) and 
  ! outputs (for time dependence)
  !===============================================================
#ifdef outGammaHI
  psys%par(:)%gammaHI = 0.0
  psys%par(:)%time = 0.0
#endif


  ! set par and src file bases for output to logfiles
  !====================================================
  fmt = "(A,'/',A,'_',I3.3)"
  write(snpbase,fmt) trim(GV%SnapPath),   trim(GV%ParFileBase),    GV%CurSnapNum
  write(srcbase,fmt) trim(GV%SourcePath), trim(GV%SourceFileBase), GV%CurSnapNum


  ! write fresh reads to the particle_data.log and source_data.log files
  !========================================================================  
  fmt = "(A,A)"
  write(str,fmt) "Fresh read from ", trim(snpbase)
  call particle_system_print_lun(psys, str, GV%pardatalun )
!  call psys%print_particle_info_lun(str,GV%pardatalun)
  write(GV%pardatalun,*)
  write(GV%pardatalun,*)
  flush(GV%pardatalun)



#ifdef useHDF5
  write(GV%srcdatalun,*) "================================================================="
  write(GV%srcdatalun,*) " HM01 G+C gammaHI for z = ", saved_gheads(GV%CurSnapNum,0)%z, ": ", &
       GV%UVB_gammaHI_cloudy
  write(GV%srcdatalun,*) "================================================================="
  write(GV%srcdatalun,*) 
#endif

  write(str,fmt) "Fresh read from ", trim(srcbase)
  call source_info_to_screen(psys,str,GV%srcdatalun)
  write(GV%srcdatalun,*)
  write(GV%srcdatalun,*)
  flush(GV%srcdatalun)
  



  ! scale the data if we need to
  !=====================================================================
  if(GV%Comoving) then
     call particle_system_scale_comoving_to_physical(psys, a, h)
!     call psys%scale_comoving_to_physical(a, h)
  endif

  ! write data after rescaling to the particle_data.log file
  !=====================================================================
  fmt = "(A,F5.3,A,F5.3,A,A)"
  write(str,fmt) "After rescaling (a=",a,",h=",h,") from ", trim(snpbase)
  call particle_system_print_lun(psys, str, GV%pardatalun ) 
!  call psys%print_particle_info_lun(str,GV%pardatalun)
  write(GV%pardatalun,*)
  write(GV%pardatalun,*)
  flush(GV%pardatalun)


  write(str,fmt) "After rescaling (a=",a,",h=",h,") from ", trim(srcbase)
  call source_info_to_screen(psys,str,GV%srcdatalun)
  write(GV%srcdatalun,*)
  write(GV%srcdatalun,*)
  flush(GV%srcdatalun)




  ! convert number density to flux for planar sources
  !==========================================================
  do i = 1,size(psys%src)
     
     if (psys%src(i)%EmisPrf == -3 .or. &
         psys%src(i)%EmisPrf == -2 .or. &
         psys%src(i)%EmisPrf == -1) then 

        ! if this is true the input luminosity is a number density 
        ! [photons/cm^3].  we want the flux that would produce this 
        ! number density in an optically thin volume
        
        write(*,*) 
        write(*,*) "  converting a photon number density to a flux"
        write(*,*) "  n_photon/cm^3                = ", psys%src(i)%L
        
        Flux = psys%src(i)%L * gconst%c * psys%box%lens_cm(1)**2
        Flux = Flux / GV%Lunit
        psys%src(i)%L = Flux
        
        write(*,*) "  photons/s from walls [1.e50] = ",  psys%src(i)%L
        write(*,*)            
        
     end if
     
  end do


  ! check test conditionals 
  !==========================================================
  if (GV%DoTestScenario) then

     if ( trim(GV%TestScenario) == "iliev_test1") then
        psys%par(:)%xHII = 1.2d-3
        psys%par(:)%xHI = 1.0d0 - psys%par(:)%xHII
        psys%par(:)%T = 1.0d4

     else if ( trim(GV%TestScenario) == "iliev_test2" ) then
        psys%par(:)%xHII = 0.0d0
        psys%par(:)%xHI = 1.0d0 - psys%par(:)%xHII
        psys%par(:)%T = 1.0d2

     else if ( trim(GV%TestScenario) == "iliev_test1He" ) then
        psys%par(:)%xHI = 1.0d0
        psys%par(:)%xHII = 0.0d0
#ifdef incHe
        psys%par(:)%xHeI = 1.0d0
        psys%par(:)%xHeII = 0.0d0
        psys%par(:)%xHeIII = 0.0d0
#endif
        psys%par(:)%T = 1.0d4

     end if
  end if



  ! set EOS particles to EOS temp if you want
  !=======================================================
#ifdef incEOS
  if (first) then
     do i = 1, size(psys%par(:))
        if (GV%EOStemp > 0.0) then
           if (psys%par(i)%eos > 0.0) then
              psys%par(i)%T = GV%EOStemp
           endif
        endif
     enddo
  endif
#endif


  ! set SFR particles to EOS temp if you want
  !=======================================================
#ifdef incSFR
  if (first) then
     do i = 1, size(psys%par(:))
        if (GV%EOStemp > 0.0) then
           if (psys%par(i)%sfr > 0.0) then
              psys%par(i)%T = GV%EOStemp
           endif
        endif
     enddo
  endif
#endif




  ! set neutral or ionized if we need to
  !=======================================================
if (first) then
   if (GV%InitxHI > 0.0) then
      psys%par(:)%xHI  = GV%InitxHI
      psys%par(:)%xHII = 1.0d0 - GV%InitxHI
      call particle_system_set_ye( psys, GV%H_mf, GV%He_mf, GV%NeBackground )
   endif
endif



  ! set constant temperature if we have one
  !=======================================================
  if (GV%IsoTemp > 0.0) psys%par(:)%T = GV%IsoTemp



  ! cap the ionization fractions and temperatures if we have to
  !================================================================
  call particle_system_enforce_x_and_T_minmax( &
       psys, GV%xfloor, GV%xceiling, GV%Tfloor, GV%Tceiling )


  ! write data after above conditionals to the particle and source log files
  !==========================================================================
  fmt = "(A,A)"
  write(str,fmt) "After test conditionals from ", trim(snpbase)
  call particle_system_print_lun(psys, str, GV%pardatalun )
  write(GV%pardatalun,*)
  write(GV%pardatalun,*)
  flush(GV%pardatalun)
 

  write(str,fmt) "After test conditionals from ", trim(srcbase)
  call source_info_to_screen(psys,str,GV%srcdatalun)
  write(GV%srcdatalun,*)
  write(GV%srcdatalun,*)
  flush(GV%srcdatalun)
  

  ! and the rest of the stuff
  !===============================================================
  GV%dt_code = PLAN%snap(GV%CurSnapNum)%RunTime / PLAN%snap(GV%CurSnapNum)%SrcRays
  call set_dt_from_dtcode( GV )
  
  GV%Tcmb_cur = gconst%t_cmb0 / a
  call get_atomic_rates(GV%Tcmb_cur, rtable, cmbT_k)

  GV%total_mass = 0.0d0
  do i = 1,size(psys%par)
     GV%total_mass = GV%total_mass + psys%par(i)%mass
  end do
  
  GV%total_lum = 0.0d0
  do i = 1,size(psys%src)
     GV%total_lum = GV%total_lum + psys%src(i)%L
  end do
  
  GV%total_atoms = GV%total_mass * GV%cgs_mass * &
       (GV%H_mf  / (gconst%protonmass) + &
       GV%He_mf / (4*gconst%protonmass) )


  GV%total_photons = (GV%TotalSimTime * GV%cgs_time / GV%LittleH) * (GV%total_lum * GV%Lunit)

  
  ! write some final data to the source log file
  !=====================================================================
  fmt = "(A,A)"
  100 format(72("-"))
  
  write(GV%srcdatalun,100)
  write(GV%srcdatalun,fmt) "Ray / Luminosity info from ", trim(srcbase)
  write(GV%srcdatalun,*) 
  write(GV%srcdatalun,'(A,ES12.5)') "dt  [code] = ", GV%dt_code
  write(GV%srcdatalun,'(A,ES12.5)') "dt  [s]    = ", GV%dt_s
  write(GV%srcdatalun,'(A,ES12.5)') "dt  [Myr]  = ", GV%dt_myr
  write(GV%srcdatalun,*)
  write(GV%srcdatalun,'(A,ES12.5)') "total photons = ", GV%total_photons
  write(GV%srcdatalun,'(A,ES12.5)') "total atoms   = ", GV%total_atoms
  write(GV%srcdatalun,'(A,ES12.5)') "photons / atoms = ", GV%total_photons / GV%total_atoms
  write(GV%srcdatalun,*)
  write(GV%srcdatalun,100)
  flush(GV%srcdatalun)
  

  call mywrite("",verb)


 

end subroutine readin_snapshot





end module main_input_mod
