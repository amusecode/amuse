module grid

  ! adaptation for amuse: public init_grid 

  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2005-04-20 (previous 2004-05-11, 2003-06-01)
  ! This module is also accepted by the F compiler (Dec 9, 2003)

  ! This module contains routines dealing with the physical 
  ! coordinate system
  ! - initcoords sets up the coordinate system
  ! - inner/outer x/y/z bound takes care of the boundary conditions

  ! Version: cartesian coordinates (x,y,z)

  ! History:
  ! 2004-05-11 - adapted to the new approach of not passing large
  !              arrays as subroutine arguments. To use the generic
  !              name vol, it is now a pointer which can point to 
  !              volx, voly, volz
  ! 2005-04-20 - adapted for 3D
  
  use file_admin, only: stdinput, log_unit, file_input
  use precision, only: dp
  use scaling, only: SCLENG
  use sizes, only: nrOfDim,neq,mbc,CART
  use my_mpi
  use mesh, only: meshx,meshy,meshz,sx,ex,sy,ey,sz,ez
  use string_manipulation, only: convert_case
  use astroconstants, only: pc,kpc,Mpc,AU

  implicit none
  private

  ! Identify type of coordinate system
  integer,parameter,public :: coordsys=CART

  ! dx,dy - cell sizes
  ! xlength,ylength - length of the entire grid
  real(kind=dp),public :: dx,dy,dz,xlength,ylength,zlength

  ! x : x-coordinate
  ! y : y-coordinate
  ! z : z-coordinate
  real(kind=dp),dimension(:),allocatable,public :: x
  real(kind=dp),dimension(:),allocatable,public :: y
  real(kind=dp),dimension(:),allocatable,public :: z

  ! volx : volume factor (for x-integration)
  ! voly : volume factor (for y-integration)
  ! volz : volume factor (for z-integration)
  ! vol  : generic name for volume, can point to volx or voly
  real(kind=dp),dimension(:,:,:),allocatable,target,public :: volx
  real(kind=dp),dimension(:,:,:),allocatable,target,public :: voly
  real(kind=dp),dimension(:,:,:),allocatable,target,public :: volz
  real(kind=dp),pointer,dimension(:,:,:),public :: vol

  real(kind=dp),dimension(2),public :: xedge
  real(kind=dp),dimension(2),public :: yedge
  real(kind=dp),dimension(2),public :: zedge

  public :: init_coords,init_grid

contains

  subroutine init_grid ( )
    
    ! Allocate the arrays
    allocate(x(sx-mbc:ex+mbc))
    allocate(y(sy-mbc:ey+mbc))
    allocate(z(sz-mbc:ez+mbc))
    ! Not used in cartesian case
    !allocate(volx(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc))
    !allocate(voly(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc))
    !allocate(volz(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc))
    
    ! point generic volume variable to volx
    vol => volx

  end subroutine init_grid

  subroutine init_coords (restart,restartfile)
    
    ! This routine initializes the coordinate variables
    
    ! This may be a fresh start or a restart of a saved run
    
    ! Case: cartesian coordinates (x,y,z)
    
    logical,intent(in) :: restart
    character(len=19),intent(in) :: restartfile

    integer :: i,j,k
    character(len=10) :: str_length_unit
    real(kind=dp) :: conversion_factor

    integer :: ierror

    if (.not.restart) then ! Fresh start
       
       ! Ask for the input if you are processor 0.
       
       if (rank == 0) then
          if (.not.file_input) write (unit=*,fmt="(a)",advance="no") &
               "2) Size of grid box (specify units): "
          read (unit=stdinput,fmt=*) xlength,ylength,zlength,str_length_unit
          write (unit=log_unit,fmt="(a,3(e10.3),a)") & 
               "2) Size of grid box : ", &
               xlength,ylength,zlength,str_length_unit
          ! Convert to cms
          call convert_case(str_length_unit,0) ! conversion to lower case
          select case (trim(adjustl(str_length_unit)))
          case ("cm","centimeter","cms","centimeters")
             conversion_factor=1.0
          case ("m","meter","ms","meters")
             conversion_factor=100.0
          case ("km","kilometer","kms","kilometers","clicks")
             conversion_factor=1000.0
          case ("pc","parsec","parsecs")
             conversion_factor=pc
          case ("kpc","kiloparsec","kiloparsecs")
             conversion_factor=kpc
          case ("mpc","megaparsec","megaparsecs")
             conversion_factor=Mpc
          case default
             write(*,*) "Length unit not recognized, assuming cm"
             conversion_factor=1.0
          end select
          xlength=xlength*conversion_factor
          ylength=ylength*conversion_factor
          zlength=zlength*conversion_factor
       endif
    else
       ! Ask for the input if you are processor 0.
       if (rank == 0) then
          call restart_grid(restartfile,xlength,ylength,zlength,ierror)
       endif
    endif

#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(xlength,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         ierror)
    call MPI_BCAST(ylength,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         ierror)
    call MPI_BCAST(zlength,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
         ierror)
#endif    

    ! Setup the coordinate grid by allocating the coordinate arrays
    call init_grid ()

    ! Scale the physical grid lengths
    xlength=xlength/scleng
    ylength=ylength/scleng
    zlength=zlength/scleng
    
    ! Cell sizes
    dx=xlength/real(max(1,meshx),dp)
    dy=ylength/real(max(1,meshy),dp)
    dz=zlength/real(max(1,meshz),dp)

    ! Fill arrays
    do i=sx-mbc,ex+mbc            
       x(i)=dx*(real(i-1,dp)+0.5_dp)
    enddo
    do j=sy-mbc,ey+mbc
       y(j)=dy*(real(j-1,dp)+0.5_dp)
    enddo
    do k=sz-mbc,ez+mbc
       z(k)=dz*(real(k-1,dp)+0.5_dp)
    enddo

    ! For higher accuracy define separate volumes for x and y
    ! This is not used in the cylindrical version
    
    !do k=sz-mbc,ez+mbc
    !   do j=sy-mbc,ey+mbc
    !      do i=sx-mbc,ex+mbc
    !         volx(i,j,k)=1.0_dp
    !         voly(i,j,k)=1.0_dp
    !         volz(i,j,k)=1.0_dp
    !      enddo
    !   enddo
    !enddo

    ! These are edge factors which are fed to the solver
    ! in [xyz]integrate. If the boundary is a singularity,
    ! it is good to set these to zero, otherwise they
    ! should be one.
    ! edge(1) = inner boundary
    ! edge(2) = outer boundary
    xedge(1:2)=1.0
    yedge(1:2)=1.0
    zedge(1:2)=1.0

  end subroutine init_coords

  !========================================================================
  subroutine restart_grid(filename,xgrid,ygrid,zgrid,ierror)
    
    ! This routine constructs the physical size of the grid
    ! (xgrid,ygrid,zgrid) from the ah3 file filename.
    ! Should be called from module coords

    use atomic, only: gamma

    character(len=19),intent(in) :: filename ! name of ah3 file
    real(kind=dp),intent(out)    :: xgrid,ygrid,zgrid ! 3D size of grid
    integer,intent(out) :: ierror

    ! AH3D header variables
    character(len=80) :: banner
    integer :: nrOfDim_in ! corresponds to parameter nrOfDim (no. of dimensions)
    integer :: neq_in     ! corresponds to parameter neq (no. of equations)
    integer :: npr_in     ! corresponds to parameter npr (no. of processors)
    integer :: refinementFactor ! not used
    integer :: nframe           ! output counter
    real(kind=dp) :: gamma_in  ! corresponds to parameter gamma (adiab. index)
    real(kind=dp) :: time      ! output time

    ! AH3D grid variables
    integer :: igrid ! counters
    integer :: xmesh,ymesh,zmesh 
    real(kind=dp) :: x_corner,y_corner,z_corner
    real(kind=dp) :: dx_in,dy_in,dz_in
    integer :: level

    ierror=0
    
    ! Read in header
    if (rank.eq.0) then
       open(unit=40,file=filename,form="UNFORMATTED",status="old")
       read(40) banner
       read(40) nrOfDim_in
       read(40) neq_in
       read(40) npr_in
       read(40) refinementFactor
       read(40) nframe
       read(40) gamma_in
       read(40) time
       
       ! Check for consistency
       if (nrOfDim_in /= nrOfDim .or. neq_in /= neq .or. npr_in /= npr .or. &
            gamma_in /= gamma ) then
          ierror=1
          write(*,*) "Error: ah3 file inconsistent with program parameters"
       endif
       
       if (ierror == 0) then
          ! Read in grids
          ! (each processor has its grid, we read in all and find the
          !  largest value to obtain the physical size of the full grid).
          xgrid=0.0
          ygrid=0.0
          zgrid=0.0
          do igrid=1,npr_in
             read(40) xmesh,ymesh,zmesh
             read(40) x_corner,y_corner,z_corner
             read(40) dx_in,dy_in,dz_in
             read(40) level
             xgrid=max(xgrid,x_corner+dx_in*(real(xmesh)-0.5))
             ygrid=max(ygrid,y_corner+dy_in*(real(ymesh)-0.5))
             zgrid=max(zgrid,z_corner+dz_in*(real(zmesh)-0.5))
          enddo
       endif

       close(40)

    endif
    
  end subroutine restart_grid

end module grid
