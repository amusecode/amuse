!>
!! \Brief This module contains data and routines for MPI parallelization
!!
!! Module for C2Ray / Capreole (3D)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2008-06-01
!!
!! \b Version: True MPI (no dummy). Also reports on OpenMP parallelization.
!! Log files for nodes 1 and higher are called 'log.n', so node 0 it is
!! 'Capreole.log'.
!! This module is also accepted by the F compiler (Dec 9, 2003)\n

module my_mpi

  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2003-06-01

  ! This module contains the routines for using the MPI library
  ! mpi_setup :   set up MPI
  ! mpi_basic :   basic MPI initialization
  ! mpi_topology: domain decomposition
  ! mpi_end:      close down the MPI interface
  ! fnd3dnbrs:    find neighbours in 3D domain decomposition
  

  ! rank 0 has a log file called Capreole.log associated with it. If
  ! the log unit is equal to 6, no file is opened and all log output
  ! is sent to standard output.
  ! The other MPI processes get their own log files, called log.1, log.2, etc.
  ! All these files are opened in mpi_setup and closed in mpi_end.
  ! Note: for very large numbers of MPI processes this is cumbersome. The
  ! log I/O should be changed for that.

  ! This is the system module:
  !use mpi

#ifdef XLF
  USE XLFUTILITY, only: hostnm => hostnm_ , flush => flush_
#endif

#ifdef IFORT
  USE IFPORT, only: hostnm, flush
#endif
  use file_admin, only: log_unit

  implicit none

  include 'mpif.h'

  integer,parameter,public :: NPDIM=3 !< dimension of problem

  integer,public :: rank            !< rank of the processor
  integer,public :: npr             !< number of processors
  integer,public :: MPI_COMM_NEW    !< the (new) communicator

  integer,dimension(NPDIM),public :: dims !< number of processors in each dimension
  integer,dimension(NPDIM),public :: grid_struct !< coordinates of the processors in the grid
  
  integer,public ::  nbrleft,nbrright  !< left and right neighbours
  integer,public ::  nbrdown,nbrup     !< up and down neighbours 
  integer,public ::  nbrabove,nbrbelow !< above and below neighbours 

  public :: mpi_setup,mpi_end, mpi_basic
  private :: mpi_topology,fnd3dnbrs

contains

  !----------------------------------------------------------------------------
  subroutine mpi_setup

    character(len=10) :: filename        ! name of the log file
    character(len=4) :: number
    integer :: ierror
    integer :: hostnm
    character(len=100) :: hostname

    ! Open processor dependent log file
    write(number,"(I4)") rank
    filename=trim(adjustl("log."//trim(adjustl(number))))
    open(unit=log_unit,file=filename,status="unknown")

    write(log_unit,*) "Log file for rank ",rank
    ! Figure out hostname
    ! NOTE: compiler dependent!!!
    ierror=hostnm(hostname)
    write(log_unit,*) "The Processor is ",trim(adjustl(hostname))
    call flush(log_unit)

    call mpi_topology

  end subroutine mpi_setup

  !----------------------------------------------------------------------------

  subroutine mpi_basic

    ! Basic MPI initialization

    integer :: mpi_ierror          ! control variable for MPI

    call MPI_INIT (mpi_ierror)  ! Initialize MPI

    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,mpi_ierror) ! Find processor rank

    ! Find total number of processors (npr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,npr,mpi_ierror)

  end subroutine mpi_basic

  !----------------------------------------------------------------------------

  subroutine mpi_topology

    ! Construct a new MPI communicator for domain decomposition

    logical,dimension(NPDIM) :: periods ! for periodic grid
    logical                  :: reorder ! reorder the MPI_COMM_WORLD
    integer          :: mpi_ierror=0

    ! Make a new topology
    dims(:)=0
    ! For plane parallel radiative transport:
    !dims(1)=1

    call MPI_Dims_create(npr,NPDIM,dims,mpi_ierror)

    periods(:)=.FALSE.      ! non-periodic boundaries

    reorder=.TRUE.
    ! makes MPI_COMM_NEW    
    call MPI_Cart_create(MPI_COMM_WORLD,NPDIM,dims,periods,reorder, &
         MPI_COMM_NEW,mpi_ierror)
    ! makes grid_struct               
    call MPI_Cart_get(MPI_COMM_NEW,NPDIM,dims, & ! makes grid_struct
         periods,grid_struct,mpi_ierror)
      
    ! Find the neighbours.
    ! My neighbors are now +/- 1 with my rank. Handle the case of the 
    ! boundaries by using MPI_PROC_NULL.
    call fnd3dnbrs ()

  end subroutine mpi_topology

  !----------------------------------------------------------------------------

  subroutine mpi_end

    ! Close MPI

    integer :: mpi_ierror=0

    ! Close log file
    close(log_unit)

    ! Close MPI
    call MPI_FINALIZE(mpi_ierror)

  end subroutine mpi_end

  !----------------------------------------------------------------------------

  subroutine fnd3dnbrs
    
    ! This routine determines the neighbours in a 3-d decomposition of
    ! the grid. This assumes that MPI_Cart_create has already been called 

    integer :: mpi_ierror=0

    call MPI_Cart_shift( MPI_COMM_NEW, 0,  1, nbrleft,  nbrright, mpi_ierror )
    call MPI_Cart_shift( MPI_COMM_NEW, 1,  1, nbrdown,  nbrup,    mpi_ierror )
    call MPI_Cart_shift( MPI_COMM_NEW, 2,  1, nbrbelow, nbrabove, mpi_ierror )

  end subroutine fnd3dnbrs

end module my_mpi
