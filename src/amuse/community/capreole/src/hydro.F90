module hydro

  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2005-04-20 (previous 2004-05-10)
  ! This module needs to be checked for the F compiler

  ! This module contains the routine which allocates the hydro
  ! arrays. This way they do not end up on stack
  
  use precision, only: dp
  use sizes, only: neq,neuler,mbc,nrOfDim,RHO,RHVX,RHVY,RHVZ,EN
  use mesh, only: sx,ex,sy,ey,sz,ez

  implicit none
  private ! variables are private by default

  real(kind=dp),dimension(:,:,:,:),allocatable,target,public :: state1
  real(kind=dp),dimension(:,:,:,:),allocatable,target,public :: state2
  real(kind=dp),dimension(:,:,:),allocatable,public   :: pressr
  real(kind=dp),dimension(:,:,:,:),allocatable,public   :: gforce

  ! These pointers are either used as a generic name (state)
  ! in some routines
  ! or as temporary storage during initialization.
  real(kind=dp),pointer,dimension(:,:,:,:),public :: state
  real(kind=dp),pointer,dimension(:,:,:,:),public :: tmpstate
  real(kind=dp),pointer,dimension(:,:,:,:),public :: stold
  real(kind=dp),pointer,dimension(:,:,:,:),public :: stnew

  integer,parameter,public :: NEW=1
  integer,parameter,public :: OLD=0

  public :: init_hydro,set_state_pointer, restart_state

contains

  !=====================================================================

  subroutine init_hydro ( )
    
    ! This routine allocates the hydrodynamic variables
    
    ! Allocate the arrays
    allocate(state1(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc,neq))
    allocate(state2(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc,neq))
    allocate(pressr(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc))
    allocate(gforce(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc,nrOfDim))
    
    ! point stnew to state2
    stold => state1
    ! point stnew to state2
    stnew => state2

    ! point generic state variable to stold
    state => stold

  end subroutine init_hydro

  !=====================================================================

  function set_state_pointer (whichone) result(state_pointer_result)

    ! Sets state pointer to either old or new state

    integer,intent(in) :: whichone
    real(kind=dp),pointer,dimension(:,:,:,:) :: state_pointer_result

    ! Point state to appropriate array
    if (whichone == NEW) then
       state_pointer_result => stnew
    else
       state_pointer_result => stold
    endif

  end function set_state_pointer

  !========================================================================

  subroutine restart_state(filename,ierror)
    
    ! This routine reads in the state variables and the time
    ! from the ah3 file filename.
    ! It should be called from the problem module.
    ! It handles multi-processor input by reading it from rank 0
    ! and sending it on the appropriate processor.
    
    use my_mpi
    use atomic, only: gamma
    
    character(len=19),intent(in) :: filename ! name of output file
    integer,intent(out) :: ierror

    ! AH3D header variables
    character(len=80) :: banner
    integer :: nrOfDim_in,neq_in,npr_in
    integer :: refinementFactor
    integer :: nframe
    real(kind=dp) :: gamma_in,time_in

    ! AH3D grid variables
    integer :: igrid ! counters
    integer :: ieq
    integer :: xmesh,ymesh,zmesh
    real(kind=dp) :: x_corner,y_corner,z_corner
    real(kind=dp) :: dx_in,dy_in,dz_in
    integer :: level
    real(kind=dp),dimension(:,:,:),allocatable :: temparray
    real(kind=dp),dimension(:,:,:,:),allocatable :: temparray2

#ifdef MPI
    integer :: status(MPI_STATUS_SIZE)
    integer :: request
    integer :: mpi_ierror
#endif

    ierror=0
    
    ! Header
    if (rank == 0) then
       open(unit=40,file=filename,form="UNFORMATTED",status="old")
       read(40) banner
       read(40) nrOfDim_in
       read(40) neq_in
       read(40) npr_in
       read(40) refinementFactor
       read(40) nframe
       read(40) gamma_in
       read(40) time_in

       ! Check for consistency
       if (nrOfDim_in /= nrOfDim .or. neq_in /= neq .or. npr_in /= npr .or. &
            gamma_in /= gamma ) then
          ierror=1
          write(30,*) "Error: ah3 file inconsistent with program parameters"
       endif
    endif

#ifdef MPI
    call MPI_BCAST(ierror,1,MPI_INTEGER,0,MPI_COMM_NEW,mpi_ierror)
#endif

    if (ierror == 0 ) then
       ! Grid: read in, but ignore. Grid should be handled by
       ! restart_grid and module coords
       if (rank == 0) then
          do igrid=1,npr_in
             read(40) xmesh,ymesh,zmesh
             read(40) x_corner,y_corner,z_corner
             read(40) dx_in,dy_in,dz_in
             read(40) level
          enddo
       endif
       
       if (rank == 0) then
          ! read in state for processor 0
          ! (this can always be done)
          read(40) state(sx:ex,sy:ey,sz:ez,RHO)
          read(40) state(sx:ex,sy:ey,sz:ez,RHVX)
          read(40) state(sx:ex,sy:ey,sz:ez,RHVY)
          read(40) state(sx:ex,sy:ey,sz:ez,RHVZ)
          read(40) state(sx:ex,sy:ey,sz:ez,EN)
          if (neq > neuler) then
             read(40) state(sx:ex,sy:ey,sz:ez,neuler+1:neq)
          endif
       endif
       
#ifdef MPI
       ! Allocate temporary arrays for reading the state variables
       allocate(temparray(ex-sx+1,ey-sy+1,ez-sz+1))
       allocate(temparray2(ex-sx+1,ey-sy+1,ez-sz+1,neuler+1:neq))
       if (rank == 0) then
          ! read in for other processors
          do igrid=2,npr_in
             do ieq=RHO,EN
                ! read in state variables one by one and send them to
                ! appropriate processor
                read(40) temparray
                call MPI_ISSEND(temparray,(ex-sx+1)*(ey-sy+1)*(ez-sz+1), &
                     MPI_DOUBLE_PRECISION,igrid-1,ieq,MPI_COMM_NEW, &
                     request,mpi_ierror)
             enddo
          enddo
          ! The non-euler variables are written in one big array, read it
          ! in and send it on
          read(40) temparray2
          call MPI_ISSEND(temparray,(ex-sx+1)*(ey-sy+1)*(ez-sz+1)*(neq-neuler), &
               MPI_DOUBLE_PRECISION,igrid-1,neuler+1,MPI_COMM_NEW, &
               request,mpi_ierror)
       else
          ! if you are a non rank 0 processor, wait for the arrays to arrive
          ! and attach it to the appropriate state variable
          do ieq=RHO,EN
             call MPI_RECV(temparray,(ex-sx+1)*(ey-sy+1)*(ez-sz+1), &
                  MPI_DOUBLE_PRECISION,0,ieq,MPI_COMM_NEW, &
                  status,mpi_ierror)
             state(sx:ex,sy:ey,sz:ez,ieq)=temparray
          enddo
          call MPI_RECV(temparray2,(ex-sx+1)*(ey-sy+1)*(ez-sz+1)*(neq-neuler), &
               MPI_DOUBLE_PRECISION,0,neuler+1,MPI_COMM_NEW, &
               status,mpi_ierror)
          state(sx:ex,sy:ey,sz:ez,neuler+1:neq)=temparray2
       endif
       deallocate(temparray)
       deallocate(temparray2)
#endif

    endif

    ! Close the input file
    if (rank == 0) close(40)

  end subroutine restart_state
  
end module hydro
