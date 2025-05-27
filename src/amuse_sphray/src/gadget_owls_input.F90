!> \file gadget_owls_input.F90

!> \brief Handles readin of GADGET OWLS/GIMIC HDF5 formatted files
!<

module gadget_owls_input_mod
use myf03_mod
use gadget_general_class
use gadget_public_header_class
use gadget_owls_header_class
use gadget_sphray_header_class
use particle_system_mod
use ion_table_class

use atomic_rates_mod, only: calc_colion_eq_fits
use global_mod, only: psys, PLAN, GV
use global_mod, only: saved_gheads

#ifdef useHDF5
use hdf5_wrapper
#endif

implicit none
private


public :: get_planning_data_gadget_owls
public :: read_Gowls_particles


contains

#ifndef useHDF5

! these are dummy subroutines so that the calls outside this file
! dont have to be wrapped with pre processor macros.

subroutine get_planning_data_gadget_owls()
  logical :: crash = .true.
  call myerr("this routine shuold not have been called","hdf5dummy",crash)
end subroutine get_planning_data_gadget_owls

subroutine read_Gowls_particles()
  logical :: crash = .true.
  call myerr("this routine shuold not have been called","hdf5dummy",crash)
end subroutine read_Gowls_particles



#else



!>   gets run planning data from Gadget OWLS/GIMIC Headers
!============================================================
subroutine get_planning_data_gadget_owls()

  character(clen), parameter :: myname = 'get_planning_data_gadget_owls_hdf5'
  logical, parameter :: crash = .true.
  integer, parameter :: verb = 2

  type(gadget_owls_header_type) :: ghead
  type(gadget_units_type) :: gunits
  type(gadget_constants_type) :: gconst

  integer(i4b) :: iSnap, fSnap    ! initial and final snapshot numbers
  integer(i4b) :: pfiles          ! files/snap for particles    
  integer(i4b) :: i,j             ! counters
  character(clen) :: snapfile     ! snapshot file name

  integer(i4b) :: loglun
  character(clen) :: logfile
 
  ! open up the planning data log file
  !======================================================
  logfile = trim(GV%OutputDir) // "/" // "particle_headers.log"
  call open_formatted_file_w(logfile,loglun)

  
  ! these global variables are read from the config file
  !======================================================
  iSnap  = GV%StartSnapNum
  fSnap  = GV%EndSnapNum    
  pfiles = GV%ParFilesPerSnap

  if ( allocated(saved_gheads) ) deallocate(saved_gheads)
  allocate( saved_gheads(iSnap:fSnap, 0:pfiles-1) )

  ! set global units and read constants
  !===================================================
  call form_gadget_snapshot_file_name(GV%SnapPath,GV%ParFileBase,iSnap,0,snapfile,hdf5bool=.true.)
  call gadget_constants_read_file( gconst, snapfile )
  call gadget_units_read_file( gunits, snapfile )
!  call gconst%read_Gowls_constants_file(snapfile)
!  call gunits%read_Gowls_units_file(snapfile)
  GV%cgs_len  = gunits%cgs_length
  GV%cgs_mass = gunits%cgs_mass
  GV%cgs_vel  = gunits%cgs_velocity
  GV%cgs_time = gunits%cgs_time
  GV%cgs_rho  = gunits%cgs_density
  GV%cgs_prs  = gunits%cgs_pressure
  GV%cgs_enrg = gunits%cgs_energy

  ! read all particle headers and write to log file
  !===================================================
  write(loglun,'(A)') "reading all GADGET OWLS/GIMIC particle header(s) ... "
  do i = iSnap,fSnap
     do j = 0,pfiles-1

        call form_gadget_snapshot_file_name(GV%SnapPath,GV%ParFileBase,i,j,snapfile,hdf5bool=.true.)
        write(loglun,'(I3,"  ",A)') i,trim(snapfile)

        call gadget_owls_header_read_file( ghead, snapfile )
        call gadget_owls_header_print_lun( ghead, loglun )
        call gadget_sphray_header_copy_owls( saved_gheads(i,j), ghead )

!        call ghead%read_Gowls_header_file(snapfile)
!        call ghead%print_Gowls_header_lun(loglun)
!        call saved_gheads(i,j)%copy_Gowls_header(ghead)

        ! make sure there is gas in this snapshot
        if (.not. ghead%npar_all(0) > 0) then
           write(*,*) "Gadget snapshot does not contain any gas particles"
           write(*,*) "Sphray cannot read dark matter particles directly, "
           write(*,*) "please calculate smoothing lengths for these particles"
           write(*,*) "and write them as gas particles. "
           stop
        end if
        
        PLAN%snap(i)%ScalefacAt = ghead%a
        PLAN%snap(i)%TimeAt = ghead%time_gyr * (gconst%SEC_PER_MEGAYEAR * 1.0d3) ! in seconds
        PLAN%snap(i)%TimeAt = PLAN%snap(i)%TimeAt * GV%LittleH / GV%cgs_time ! in code units
        
     end do
  end do

  ! close headers log file
  !========================
  close(loglun)

  ! use one header to set global variables
  !--------------------------------------------
  GV%BoxLwrs(:) = 0.0d0
  GV%BoxUprs(:) = saved_gheads(iSnap,0)%boxlen
  
  GV%OmegaM = saved_gheads(iSnap,0)%OmegaM
  GV%OmegaL = saved_gheads(iSnap,0)%OmegaL
  GV%OmegaB = saved_gheads(iSnap,0)%OmegaB 
  
  GV%LittleH = saved_gheads(iSnap,0)%h
  
  ! write units to log file
  !===================================================

  logfile = trim(GV%OutputDir) // "/" // "code_units.log"
  call open_formatted_file_w(logfile,loglun)
  call gadget_units_print_lun( gunits, loglun, saved_gheads(iSnap,0)%h ) 
  close(loglun)

end subroutine get_planning_data_gadget_owls







!> reads a Gadget OWLS/GIMIC HDF5 snapshot into a particle array  
!========================================================================
subroutine read_Gowls_particles()

  character(clen), parameter :: myname="read_Gowls_particles" 
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt
  
  real(r4b), allocatable :: rblck(:)
  real(r4b), allocatable :: rblck3(:,:)
  integer(i4b), allocatable :: iblck(:)
  integer(i8b) :: ngasread

  type(gadget_owls_header_type) :: ghead
  type(gadget_sphray_header_type) :: shead
  character(clen) :: snapfile, VarName, GroupName
  integer(i4b) :: fh
  integer(i8b) :: i
  integer(i4b) :: err

  integer(i8b) :: npar, ngas, nmass
  integer(i8b) :: npar1, ngas1, nmass1
  logical :: varmass(0:5)
  integer(i4b) :: fn

  real(r8b) :: meanweight
  logical :: caseA(2)
  real(r8b) :: MB

  real(r8b) :: Tdum
  real(r8b) :: Hmf
  real(r8b) :: nH8
  real(r8b) :: T8
  real(r8b) :: xvec(5)

  type(gadget_constants_type) :: gconst
  type(ion_table_type) :: itab
  real(r8b) :: redshift

  logical :: hdf5bool

  ! set hdf5 boolean
  !======================================================
  hdf5bool = .true.

  ! set local particle numbers
  !============================
  shead = saved_gheads( GV%CurSnapNum, 0 )
  varmass = (shead%npar_all > 0 .and. shead%mass == 0)
  npar = sum(shead%npar_all)
  ngas = shead%npar_all(0)
  nmass = sum(shead%npar_all, mask=varmass)


  ! do Gadget dummy checks
  !============================
  if (ngas == 0) call myerr("snapshot has no gas particles",myname,crash)

  ! calculate bytes per particle and allocate particle array
  !===========================================================
  MB = GV%bytesperpar * real(ngas) / 2.0d0**20
  GV%MB = GV%MB + MB

  fmt="(A,F10.4,A,I10,A)"
  write(str,fmt) "   allocating ", MB, " MB for ", ngas, " particles"
  call mywrite(str,verb) 

  allocate (psys%par(ngas), stat=err)
  if (err /= 0) call myerr("failed to allocate par",myname,crash)


  ! now read all snapshot files
  !==============================          
  ngasread = 0
  GroupName = 'PartType0/'
  files: do fn = 0, shead%nfiles-1


     ! recall the header info
     !-----------------------------------------------------------!  
     shead   = saved_gheads( GV%CurSnapNum, fn )
     varmass = (shead%npar_file > 0 .and. shead%mass == 0)
     npar1   = sum(shead%npar_file)
     ngas1   = shead%npar_file(0)
     nmass1  = sum(shead%npar_file, mask=varmass)
     if (ngas1 == 0) cycle


     ! begin read
     !-----------------------------------------------------------!  
     call form_gadget_snapshot_file_name(GV%SnapPath,GV%ParFileBase,GV%CurSnapNum,fn,snapfile,hdf5bool)
     call mywrite("   reading particle snapshot file: "//trim(snapfile), verb)
     call hdf5_open_file(fh, snapfile, readonly=.true.)
     call gadget_owls_header_read_lun(ghead,fh)



     ! read positions 
     !-----------------------------------------------------------!  
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for pos",myname,crash)
     VarName = 'Coordinates'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck3)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(1) = rblck3(1,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(2) = rblck3(2,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(3) = rblck3(3,i)
     deallocate(rblck3)
     call gadget_data_attributes_read_lun( pos_attrs, fh, GroupName, VarName )



     ! read velocities 
     !-----------------------------------------------------------!  
#ifdef incVel
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for vel",myname,crash)
     VarName = 'Velocity'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck3)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(1) = rblck3(1,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(2) = rblck3(2,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(3) = rblck3(3,i)
     deallocate(rblck3)
     call gadget_data_attributes_read_lun( vel_attrs, fh, GroupName, VarName )
!     call vel_attrs%read_lun( fh, GroupName, VarName )
#endif

     ! read id's 
     !-----------------------------------------------------------!  
     allocate(iblck(ngas1), stat=err )
     if(err/=0) call myerr("allocating iblck for ID",myname,crash)
     VarName = 'ParticleIDs'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),iblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%id = iblck(i)
     deallocate(iblck)
     call gadget_data_attributes_read_lun( id_attrs, fh, GroupName, VarName ) 
!    call id_attrs%read_lun( fh, GroupName, VarName )

     ! read masses 
     !-----------------------------------------------------------!  

     ! if gas particles are variable mass
     if (varmass(0)) then  
        allocate(rblck(ngas1), stat=err)
        if(err/=0) call myerr("allocating rblck for mass",myname,crash)
        VarName = 'Mass'
        call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
        forall(i=1:ngas1) psys%par(ngasread+i)%mass = rblck(i)
        deallocate(rblck)
        call gadget_data_attributes_read_lun( mass_attrs, fh, GroupName, VarName ) 
!       call mass_attrs%read_lun( fh, GroupName, VarName )

     ! if gas particles are isomass
     else
        psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(0)       
     end if


     ! read temperature
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for T",myname,crash)
     VarName = 'Temperature'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%T = rblck(i)
     deallocate(rblck)
     call gadget_data_attributes_read_lun( T_attrs, fh, GroupName, VarName )
!     call T_attrs%read_lun( fh, GroupName, VarName ) 

     ! read EOS
     !-----------------------------------------------------------!  
#ifdef incEOS
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for EOS",myname,crash)
     VarName = 'OnEquationOfState'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%eos = rblck(i)
     deallocate(rblck)
     call gadget_data_attributes_read_lun( eos_attrs, fh, GroupName, VarName )
!     call eos_attrs%read_lun( fh, GroupName, VarName )
#endif

     ! read density 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for rho",myname,crash)
     VarName = 'Density'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%rho = rblck(i)
     deallocate(rblck)
     call gadget_data_attributes_read_lun( rho_attrs, fh, GroupName, VarName ) 
!    call rho_attrs%read_lun( fh, GroupName, VarName )

     ! read smoothing lengths 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for hsml",myname,crash)
     VarName = 'SmoothingLength'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%hsml = rblck(i)
     deallocate(rblck)
     call gadget_data_attributes_read_lun( hsml_attrs, fh, GroupName, VarName )
!     call hsml_attrs%read_lun( fh, GroupName, VarName )

     ! read Hydrogen mass fractions
     !-----------------------------------------------------------!  
#ifdef incHmf
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for Hmf",myname,crash)
     VarName = 'ElementAbundance/Hydrogen'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%Hmf = rblck(i)
     deallocate(rblck)
     call gadget_data_attributes_read_lun( Hmf_attrs, fh, GroupName, VarName )
!     call Hmf_attrs%read_lun( fh, GroupName, VarName )
#endif


     ! read Helium mass fractions
     !-----------------------------------------------------------!  
#ifdef incHemf
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for Hemf",myname,crash)
     VarName = 'ElementAbundance/Helium'
     call hdf5_read_data(fh,trim(GroupName)//trim(VarName),rblck)
     forall(i=1:ngas1) psys%par(ngasread+i)%Hemf = rblck(i)
     deallocate(rblck)
     call gadget_data_attributes_read_lun( Hemf_attrs, fh, GroupName, VarName )
!     call Hemf_attrs%read_lun( fh, GroupName, VarName )
#endif


     ngasread = ngasread + ngas1
     call hdf5_close_file(fh)


  end do files


  ! set caseA true or false for collisional equilibrium
  !-----------------------------------------------------
  caseA = .false.
  if (GV%HydrogenCaseA) caseA(1) = .true.
  if (GV%HeliumCaseA)   caseA(2) = .true.


  ! calculate xHI from CLOUDY iontables
  !-----------------------------------------------------------!  
  call mywrite('',verb) 
  call mywrite("   calculating input xHI from CLOUDY tables",verb)

  call read_ion_table_file( "../data/ionization_tables/h1.hdf5", itab )
  redshift = ghead%z

  ! first get the gammaHI from the uniform UVB at this redshift
  GV%UVB_gammaHI_cloudy = return_gammaHI_at_z( itab, redshift )

 
  ! loop through the gas particles and interpolate from the table
  !---------------------------------------------------------------
  do i = 1,ngas

     ! set individual Hydrogen mass fractions if we have them
#ifdef incHmf
     Hmf = psys%par(i)%Hmf
#else
     Hmf = GV%H_mf
#endif


     ! get Hydrogen number density and temperature
     nH8 = psys%par(i)%rho * GV%cgs_rho * ghead%h**2 / ghead%a**3 * &
           Hmf / gconst%PROTONMASS
     T8 = psys%par(i)%T

     ! if we have EOS particles set their temperature if we need to
#ifdef incEOS
     if (psys%par(i)%eos > 0.0) then
        if (GV%EOStemp > 0.0) then
           T8 = GV%EOStemp
           psys%par(i)%T = T8
        endif
     endif
#endif

     psys%par(i)%xHI = &
          interpolate_ion_table( itab, redshift, log10(T8), log10(nH8) )
     
  end do



#ifdef incCloudy
  psys%par(:)%xHI_cloudy = psys%par(:)%xHI
  fmt = "(T7, A, 2ES15.5)"
  write(str,fmt) "min/max xHI_cloudy = ", minval( psys%par%xHI_cloudy ), maxval( psys%par%xHI_cloudy )
  call mywrite(str,verb) 
#endif
  call mywrite('',verb)

  

  ! set xHII from xHI 
  !-----------------------------------------------------------!  
  psys%par%xHII = 1.0d0 - psys%par%xHI



  ! if Helium, initialize ionization fractions to collisional equilibrium
  !------------------------------------------------------------------------
#ifdef incHe
  call psys%set_ci_eq(caseA, DoH=.false., DoHe=.true., fit="hui")
#endif


  ! set the electron fractions from the ionization fractions
  !----------------------------------------------------------
  call particle_system_set_ye( psys, GV%H_mf, GV%He_mf, GV%NeBackground )



end subroutine read_Gowls_particles




#endif


end module gadget_owls_input_mod
