module times

  ! Module for Capreole (f90)
  ! Author: Garrelt Mellema
  ! Date: 2004-05-11
  ! This module needs to be checked for the F compiler

  ! This module contains the routine which allocates the hydro
  ! arrays. This way they do not end up on stack
  
  use file_admin, only: stdinput, log_unit, file_input
  use precision, only: dp
  use my_mpi
  use scaling, only: SCTIME
  use string_manipulation, only: convert_case
  use astroconstants, only: YEAR

  implicit none
  private

  real(kind=dp),public :: frametime
  integer,public       :: LastFrame
  real(kind=dp),public :: time
  real(kind=dp),public :: dt

  public :: init_time

contains

  subroutine init_time (restart,restartfile)
    
    ! This routine initializes all time variables
    
    ! This may be a fresh start or a restart of a saved run
    
    logical,intent(in) :: restart
    character(len=19),intent(in) :: restartfile

    character(len=10) :: str_time_unit
    integer :: ierror

    if (restart) then 
       
       call restart_time(restartfile,time,ierror)
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,0,&
            MPI_COMM_NEW,ierror)
#endif

    else ! Fresh start
       
       ! Initialize the time to zero
       time=0.0d0
    endif

    ! Initialize the time step to zero
    dt=0.0d0

    ! Ask for the input if you are processor 0.
    if (rank == 0) then
       if (.not.file_input) then
          write (*,"(//,A,/)") "----- Output -----"
          write (*,"(A,$)") "1) Time between outputs (specify unit): "
       endif
       read (stdinput,*) frametime,str_time_unit
       if (.not.file_input) write (*,"(A,$)") "2) How many output frames: "
       read (stdinput,*) LastFrame
    endif
    
    ! report input parameters
    if (rank == 0) then
       write (log_unit,"(//,A,/)") "----- Output -----"
       write (log_unit,"(A,1PE10.3,A)") "1) Time between outputs: ", & 
            frametime,str_time_unit
       write (log_unit,"(A,I4)") "2) Number of output frames: ",LastFrame
       ! Convert to seconds
       call convert_case(str_time_unit,0) ! conversion to lower case
       select case (trim(adjustl(str_time_unit)))
       case ("s","sec","second","secs","seconds")
       case ("years","yrs","year","yr")
          frametime=frametime*YEAR
       case ("myears","myrs","myear","myr")
          frametime=frametime*1e6*YEAR
       case ("gyears","gyrs","gyear","gyr")
          frametime=frametime*1e9*YEAR
       case default
          write(log_unit,*) "Time unit not recognized, assuming seconds"
       end select
    endif
    
    
#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(frametime,1,MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,ierror)
    call MPI_BCAST(LastFrame,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
#endif
    
    ! Scale time
    frametime=frametime/sctime
    
  end subroutine init_time

  !========================================================================
  subroutine restart_time(filename,simtime,ierror)
    
    ! This routine retrieves the time (simtime)
    ! from the ah3 file filename.
    ! Should be called from module time

    character(len=19),intent(in) :: filename ! name of ah3 file
    real(kind=dp),intent(out)    :: simtime ! time of restart
    integer,intent(out) :: ierror

    ! AH3D header variables
    character(len=80) :: banner
    integer :: nrOfDim_in ! corresponds to parameter nrOfDim (no. of dimensions)
    integer :: neq_in     ! corresponds to parameter neq (no. of equations)
    integer :: npr_in     ! corresponds to parameter npr (no. of processors)
    integer :: refinementFactor ! not used
    integer :: nframe           ! output counter
    real(kind=dp) :: gamma_in  ! corresponds to parameter gamma (adiab. index)

    ierror=0
    
    ! Read in header
    if (rank == 0) then
       open(unit=40,file=filename,form="UNFORMATTED",status="old")

       read(40) banner
       read(40) nrOfDim_in
       read(40) neq_in
       read(40) npr_in
       read(40) refinementFactor
       read(40) nframe
       read(40) gamma_in
       read(40) simtime
       
       close(40)
    endif
    
    write(log_unit,*) "Sim time is", simtime/YEAR," years"
    
    ! Scale time
    simtime=simtime/SCTIME
    
  end subroutine restart_time

end module times
