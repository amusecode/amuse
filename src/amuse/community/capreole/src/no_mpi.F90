!>
!! \brief This module contains data and routines for MPI parallelization
!!
!! Module for C2Ray / Capreole (3D)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2003-06-01
!!
!! \b Version: This is a dummy module for systems where there is no MPI for F90.
!!
!! This module is also accepted by the F compiler (Dec 9, 2003)\n

module my_mpi

  use file_admin, only: log_unit, results_dir

#ifdef XLF
  USE XLFUTILITY, only: hostnm => hostnm_ , flush => flush_
#endif

#ifdef IFORT 
  USE IFPORT, only: hostnm, flush
#ifdef _OPENMP
  USE OMP_LIB, only: omp_get_num_threads, omp_get_thread_num
#endif
#endif

  implicit none

  integer,parameter,public :: NPDIM=3 !< dimension of problem

  ! All of these are set to be consistent with the MPI version
  integer,public :: rank              !< rank of the processor
  integer,public :: npr               !<< number of processors
  integer,public :: nthreads        !< number of threads (per processor)
  integer,public :: MPI_COMM_NEW      !< the (new) communicator (dummy)

  integer,dimension(NPDIM),public :: dims !< number of processors in each dimension (dummy)
  integer,dimension(NPDIM),public :: grid_struct !< coordinates of the processors in the grid (dummy)

  integer,public ::  nbrleft,nbrright  !< left and right neighbours
  integer,public ::  nbrdown,nbrup     !< up and down neighbours 
  integer,public ::  nbrabove,nbrbelow !< above and below neighbours 

  integer,parameter,public :: MPI_PROC_NULL=-1

#ifdef SUN
  integer :: hostnm
#endif

  public :: mpi_setup,mpi_end
  private :: mpi_basic,mpi_topology,fnd2dnbrs

contains

  !----------------------------------------------------------------------------

  subroutine mpi_setup ( )

    character(len=512) :: filename        ! name of the log file
    character(len=4) :: number
    integer :: ierror
    integer :: tn
    character(len=100) :: hostname

    call mpi_basic ()

    if (log_unit /= 6) then
       filename=trim(adjustl(trim(adjustl(results_dir))//"Capreole.log"))
       open(unit=log_unit,file=filename,status="unknown",action="write", &
            position="append")
    endif
    write(unit=log_unit,fmt="(A)") "Log file for Capreole run"

    nthreads=1
    ! Figure out hostname
    ! NOTE: compiler dependent!!!
    ierror=hostnm(hostname)
    if (ierror == 0) then
       write(log_unit,*) "Running on processor named ",hostname
    else 
       write(log_unit,*) "Error establishing identity of processor."
    endif

    ! Report number of OpenMP threads
    !$omp parallel default(shared)
#ifdef _OPENMP
    nthreads=omp_get_num_threads()
#endif
    !$omp end parallel
#ifdef _OPENMP
    write(log_unit,*) ' Number of OpenMP threads is ',nthreads
#endif

    ! Let OpenMP threads report
    !$omp parallel default(private)
#ifdef _OPENMP
    tn=omp_get_thread_num()+1
    write(log_unit,*) 'Thread number ',tn,' reporting'
#endif
    !$omp end parallel

!    flush(log_unit)

    call mpi_topology ()

  end subroutine mpi_setup

  !----------------------------------------------------------------------------

  subroutine mpi_basic ( )

    rank=0 ! Find processor rank

    ! Find total number of processors (npr)
    npr=1

  end subroutine mpi_basic

  !----------------------------------------------------------------------------

  subroutine mpi_topology ( )


    ! Make a new topology
    dims=1

    ! makes MPI_COMM_NEW    
    MPI_COMM_NEW=0

    ! makes grid_struct               
    grid_struct=0
      
    ! Find the neighbours.
    ! My neighbors are now +/- 1 with my rank. Handle the case of the 
    ! boundaries by using MPI_PROC_NULL.
    call fnd2dnbrs ( )

  end subroutine mpi_topology

  !----------------------------------------------------------------------------

  subroutine mpi_end ( )

    ! Close log file
    close(log_unit)

  end subroutine mpi_end

  !----------------------------------------------------------------------------

  subroutine fnd2dnbrs ( )
    
    ! This routine determines the neighbours in a 2-d decomposition of
    ! the grid. This assumes that MPI_Cart_create has already been called 

    ! Single processor version

    nbrleft=MPI_PROC_NULL
    nbrright=MPI_PROC_NULL
    nbrup=MPI_PROC_NULL
    nbrdown=MPI_PROC_NULL
    nbrabove=MPI_PROC_NULL
    nbrbelow=MPI_PROC_NULL
    
  end subroutine fnd2dnbrs

end module my_mpi
