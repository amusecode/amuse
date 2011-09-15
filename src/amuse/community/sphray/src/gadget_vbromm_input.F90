!> \file gadget_vbromm_input.F90

!> \brief Handles readin of GADGET 2.0 raw binary formatted files
!! from Volker Bromm's group
!<

module gadget_vbromm_input_mod
use myf03_mod
use gadget_general_class
use gadget_public_header_class
use gadget_sphray_header_class
use particle_system_mod, only: particle_system_type

use global_mod, only: psys, PLAN, GV
use global_mod, only: saved_gheads
implicit none
private

public :: get_planning_data_gadget_vbromm
public :: read_Gvbromm_particles
public :: set_temp_from_u_vbromm



contains


!>  gets run planning data from Gadget VBromm Headers
!======================================================
subroutine get_planning_data_gadget_vbromm()

  type(gadget_public_header_type) :: ghead
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

  ! set global units
  !===================================================
  GV%cgs_len  = gunits%cgs_length
  GV%cgs_mass = gunits%cgs_mass
  GV%cgs_vel  = gunits%cgs_velocity
  GV%cgs_time = gunits%cgs_time
  GV%cgs_rho  = gunits%cgs_density
  GV%cgs_prs  = gunits%cgs_pressure
  GV%cgs_enrg = gunits%cgs_energy  

  ! read all particle headers and write to log file
  !===================================================
  write(loglun,'(A)') "reading all GADGET Vbromm particle header(s) ... "
  do i = iSnap,fSnap
     do j = 0,pfiles-1

        call form_gadget_snapshot_file_name(GV%SnapPath, GV%ParFileBase, i, j, snapfile, hdf5bool=.false.)
        write(loglun,'(I3,"  ",A)') i, trim(snapfile)

        call gadget_public_header_read_file( ghead, snapfile )
        call gadget_public_header_print_lun( ghead, loglun )
        call gadget_sphray_header_copy_public( saved_gheads(i,j), ghead )

!        call ghead%read_Gpublic_header_file(snapfile)
!        call ghead%print_Gpublic_header_lun(loglun)
!        call saved_gheads(i,j)%copy_Gpublic_header(ghead)

        saved_gheads(i,j)%OmegaB = 0.045 ! this should be moved to the config file
                                         ! but its not used in the code now

        saved_gheads(i,j)%time_gyr = gadget_public_header_return_gyr(ghead) 
!        saved_gheads(i,j)%time_gyr = ghead%return_gyr()

        ! make sure there is gas in this snapshot
        !-------------------------------------------
        if (.not. ghead%npar_all(0) > 0) then
           write(*,*) "Gadget snapshot does not contain any gas particles"
           write(*,*) "Sphray cannot read dark matter particles directly, "
           write(*,*) "please calculate smoothing lengths for these particles"
           write(*,*) "and write them as gas particles. "
           stop
        end if

        if (GV%Comoving) then
           PLAN%snap(i)%ScalefacAt = ghead%a
           PLAN%snap(i)%TimeAt = gadget_public_header_return_gyr(ghead) * (gconst%SEC_PER_MEGAYEAR * 1.0d3) ! in seconds
           PLAN%snap(i)%TimeAt = PLAN%snap(i)%TimeAt * GV%LittleH / GV%cgs_time ! in code units
        else
           PLAN%snap(i)%ScalefacAt = 1.0d0 / (1.0d0 + ghead%z)
           PLAN%snap(i)%TimeAt = ghead%a
        end if

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

end subroutine get_planning_data_gadget_vbromm


!> reads a Gadget snapshot from Volker Bromm's group
!=================================================================
subroutine read_Gvbromm_particles()

  character(clen), parameter :: myname="read_Gvbromm_particles"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt 

  real(r4b), allocatable :: rblck(:)
  real(r4b), allocatable :: rblck3(:,:)
  integer(i4b), allocatable :: iblck(:)
  integer(i8b) :: ngasread

  type(gadget_public_header_type) :: ghead
  type(gadget_sphray_header_type) :: shead
  character(clen) :: snapfile
  integer(i4b) :: lun,i
  integer(i4b) :: err

  integer(i8b) :: npar, ngas, nmass
  integer(i8b) :: npar1, ngas1, nmass1
  logical :: varmass(0:5)
  integer(i4b) :: fn

  real(r8b) :: meanweight
  logical :: caseA(2)
  real(r8b) :: xvec(5)
  real(r8b) :: Tdum
  real(r8b) :: MB 
  real(r8b) :: nH_over_nHe

  logical :: hdf5bool

  ! set hdf5 boolean
  !======================================================
  hdf5bool = .false.

  ! set local particle numbers
  !============================
  shead = saved_gheads( GV%CurSnapNum, 0 )
  varmass = (shead%npar_all > 0 .and. shead%mass == 0)
  npar = sum(shead%npar_all)
  ngas = shead%npar_all(0)
  nmass = sum(shead%npar_all, mask=varmass)

  ! do Gadget dummy checks
  !============================
  if (ngas .EQ. 0) call myerr("snapshot has no gas particles",myname,crash)

  ! calculate bytes per particle and allocate particle array
  !===========================================================
  MB = GV%bytesperpar * real(ngas) / 2.0d0**20
  GV%MB = GV%MB + MB

  fmt="(A,F10.4,A,I10,A)"
  write(str,fmt) "allocating ", MB, " MB for ", ngas, " particles"
  call mywrite(str,verb)

  allocate (psys%par(ngas), stat=err)
  if (err /= 0) call myerr("failed to allocate par",myname,crash)



  ! now read all snapshot files
  !==============================          
  ngasread = 0
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
     call mywrite("reading vbromm gadget snapshot file "//trim(snapfile), verb)
     call open_unformatted_file_r( snapfile, lun )
     call gadget_public_header_read_lun( ghead, lun )



     ! read positions 
     !-----------------------------------------------------------!  
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for pos",myname,crash)
     read(lun, iostat=err) rblck3
     if (err/=0) call myerr("reading rblk3 for pos",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(1) = rblck3(1,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(2) = rblck3(2,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%pos(3) = rblck3(3,i)
     deallocate(rblck3)

     ! read velocities 
     !-----------------------------------------------------------!  
#ifdef incVel
     allocate(rblck3(3,ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck3 for vel",myname,crash)
     read(lun, iostat=err) rblck3
     if (err/=0) call myerr("reading rblk3 for vel",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(1) = rblck3(1,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(2) = rblck3(2,i)
     forall(i=1:ngas1) psys%par(ngasread+i)%vel(3) = rblck3(3,i)
     deallocate(rblck3)
#else
     read(lun, iostat=err)
     if (err/=0) call myerr("dummy reading vel",myname,crash) 
#endif


     ! read id's 
     !-----------------------------------------------------------!  
     allocate(iblck(ngas1), stat=err )
     if(err/=0) call myerr("allocating iblck for ID",myname,crash)
     read(lun, iostat=err) iblck  
     if (err/=0) call myerr("reading iblk for ID",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%id = iblck(i)
     deallocate(iblck)

     ! read masses 
     !-----------------------------------------------------------!  

     ! if there are variable mass particles of any type
     if (nmass1 > 0) then  

        ! if there are variable mass gas particles
        if (varmass(1)) then  
           allocate(rblck(ngas1), stat=err)
           if(err/=0) call myerr("allocating rblck for mass",myname,crash)
           read(lun, iostat=err) rblck 
           if (err/=0) call myerr("reading rblk for mass",myname,crash) 
           forall(i=1:ngas1) psys%par(ngasread+i)%mass = rblck(i)
           deallocate(rblck)

        ! just dummy read non gas particles
        else 
           read(lun)  
           psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(1)
        end if

     ! if none of the particles are variable mass
     else  
        psys%par(ngasread+1:ngasread+ngas1)%mass = ghead%mass(1)
     end if


     ! read temperature (internal energy / unit mass for now)
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for u",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for u",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%T = rblck(i) 
     deallocate(rblck)


     ! read density 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for rho",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for rho",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%rho = rblck(i)
     deallocate(rblck)


     ! read smoothing lengths 
     !-----------------------------------------------------------!  
     allocate(rblck(ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for hsml",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for hsml",myname,crash) 
     forall(i=1:ngas1) psys%par(ngasread+i)%hsml = rblck(i)
     deallocate(rblck)

     ! read ionization fractions
     !-----------------------------------------------------------!  
     
     ! CAFG: First, big array containing all the abundances.
     allocate(rblck(10*ngas1), stat=err)
     if(err/=0) call myerr("allocating rblck for all abundances",myname,crash)
     read(lun, iostat=err) rblck
     if (err/=0) call myerr("reading rblk for all abundances",myname,crash) 
     
     ! CAFG: Here is where we should assign the desired abundances to the
     ! relevant particle data.
     forall(i=1:ngas1) psys%par(ngasread+i)%xHII = rblck(10*(i-1)+2)
     
 
#ifdef incHe
     ! read HeII = nHeII / nH
     !-----------------------------------------------------------!  
     forall(i=1:ngas1) psys%par(ngasread+i)%xHeII = rblck(10*(i-1)+9)


     ! read HeIII = nHeIII / nH
     !-----------------------------------------------------------!  
     forall(i=1:ngas1) psys%par(ngasread+i)%xHeIII = rblck(10*(i-1)+10)

#endif
     
     deallocate(rblck)

     ngasread = ngasread + ngas1
     close(lun)

  end do files



  ! set xHI from xHII 
  !----------------------------------------------------------------!  
  psys%par(:)%xHI = 1.0d0 - psys%par(:)%xHII


  ! use Hydrogen and Helium mass fractions from config file to
  ! set the Helium ionization fractions
  !----------------------------------------------------------------!  
#ifdef incHe
  nH_over_nHe = 4 * GV%H_mf / GV%He_mf
  psys%par(:)%xHeII  = psys%par(:)%xHeII  * nH_over_nHe
  psys%par(:)%xHeIII = psys%par(:)%xHeIII * nH_over_nHe
  psys%par(:)%xHeI   = 1.0d0 - psys%par(:)%xHeII - psys%par(:)%xHeIII
#endif 

  ! set ye from ionization fractions
  !----------------------------------------------------------------!  
  psys%par(:)%ye = psys%par(:)%xHII
#ifdef incHe
  psys%par(:)%ye = psys%par(:)%ye + psys%par(:)%xHeII + 2 * psys%par(:)%xHeIII
#endif

  ! convert the internal energy to temperature using ye 
  !----------------------------------------------------------------!  
  
  call set_temp_from_u_vbromm(psys, GV%H_mf, GV%cgs_enrg, GV%cgs_mass)


  ! shrink particles with negative IDs to 
  !----------------------------------------------------------------!  
  where( psys%par(:)%id < 0 ) psys%par(:)%hsml = 0.0




end subroutine read_Gvbromm_particles


!> converts internal energies / unit mass to temperature K using Hmf and ye = ne/nH
!=======================================================================================
subroutine set_temp_from_u_vbromm(psys, dfltH_mf, cgs_enrg, cgs_mass)

  type(particle_system_type) :: psys
  real(r8b), intent(in) :: dfltH_mf
  real(r8b), intent(in) :: cgs_enrg
  real(r8b), intent(in) :: cgs_mass
  integer(i8b) :: i
  real(r8b) :: Hmf
  real(r8b) :: mu
  real(r8b) :: Tdum
  type(gadget_constants_type) :: gconst

  do i = 1,size(psys%par)

#ifdef incHmf
     Hmf = psys%par(i)%Hmf
#else
     
     Hmf = dfltH_mf    

#endif

     mu = 4.0d0 / (3.0d0 * Hmf + 1.0d0 + 4.0d0 * Hmf * psys%par(i)%ye)

     ! CAFG: artificially 'correct' abnormally large values, which are
     ! producing floating overflow
     if (psys%par(i)%T .GT. 1000000000.0d0) psys%par(i)%T = 1000.0d0

     Tdum = mu * gconst%PROTONMASS / gconst%BOLTZMANN * (gconst%GAMMA - 1.0d0) * psys%par(i)%T
     psys%par(i)%T = Tdum * cgs_enrg / cgs_mass

  end do


end subroutine set_temp_from_u_vbromm


end module gadget_vbromm_input_mod
