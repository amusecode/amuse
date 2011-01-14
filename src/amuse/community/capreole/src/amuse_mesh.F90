module mesh
  
  ! adaptation for amuse (public fnd3ddecomp)

  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2005-04-20 (prev 2003-06-01)
  ! This module is also accepted by the F compiler (Dec 9, 2003)
  !
  ! This module contains the routines to set up a computational mesh,
  ! which can be distributed over several processors using
  ! MPI routines.
  !
  ! 3D version

  use precision, only: dp
  use my_mpi
  use file_admin, only: stdinput,log_unit,ah3,file_input

  implicit none

  private

  integer,public :: meshx !< total mesh size in x, y, z
  integer,public :: meshy !< total mesh size in x, y, z
  integer,public :: meshz !< total mesh size in x, y, z
  integer,public :: sx !< local mesh start and end coordinates
  integer,public :: ex !< local mesh start and end coordinates
  integer,public :: sy !< local mesh start and end coordinates
  integer,public :: ey !< local mesh start and end coordinates
  integer,public :: sz !< local mesh start and end coordinates
  integer,public :: ez !< local mesh start and end coordinates

  public :: init_mesh !< subroutine to initialize mesh
  public :: fnd3ddecomp

#ifdef MPI
  integer :: mpi_ierror !< mpi error flag
#endif

contains

  !----------------------------------------------------------------------------

  !> Read in the mesh dimensions, and distribute them over the
  !! processors
  subroutine init_mesh (restart,restartfile)

    logical,intent(in) :: restart !< restart or not
    character(len=*),intent(in) :: restartfile !< name of restart file

    integer :: ierror=0 !< error flag

    ! Ask for the input if you are processor 0.
    if (rank == 0) then
       if (.not.restart) then ! Fresh start
          
          ! Ask for input
          if (.not.file_input) then
             print "(2/,A,/)", "----- Mesh -----"
             write(unit=*,fmt="(a)",advance="no") "1) Number of mesh points: "
          endif
          read (unit=stdinput,fmt=*) meshx,meshy,meshz
          
          ! Report back mesh
          write(unit=log_unit,fmt="(2/,A,/)") "----- Mesh -----"
          write(unit=log_unit,fmt="(A,3I5)") "1) Number of mesh points: ", &
               meshx,meshy,meshz
       else
          
          ! Set mesh from restart file
          call restart_mesh(restartfile,meshx,meshy,meshz,ierror)
          
       endif

    endif

#ifdef MPI
    ! Distribute the total mesh size over all processors
    call MPI_BCAST(meshx,1,MPI_INTEGER,0,MPI_COMM_NEW,mpi_ierror)
    call MPI_BCAST(meshy,1,MPI_INTEGER,0,MPI_COMM_NEW,mpi_ierror)
    call MPI_BCAST(meshz,1,MPI_INTEGER,0,MPI_COMM_NEW,mpi_ierror)
#endif

    ! Compute the decomposition of the mesh over the processors
    ! (find 3D decomposition)
    ! sx = start x coordinate, ex = end x coordinate
    ! sy = start y coordinate, ey = end y coordinate
    ! sz = start z coordinate, ez = end z coordinate
    call fnd3ddecomp ()

    ! Report the mesh for the local processor
    write(unit=log_unit,fmt=*) "Local mesh: ",sx,ex,sy,ey,sz,ez
    
  end subroutine init_mesh

  !> This routine retrieved the mesh size
  !! (xmesh,ymesh,zmesh) from the ah3 file 'filename'.
  !! Should be called from module mesh
  subroutine restart_mesh (filename,xmesh,ymesh,zmesh,ierror)
    
    use sizes, only: nrOfDim, neq
    use atomic, only: gamma

    character(len=*),intent(in) :: filename !< name of ah3 file
    integer,intent(out) :: xmesh,ymesh,zmesh !< 3D size of mesh
    integer,intent(out) :: ierror !< error flag

    ! AH3D header variables
    character(len=80) :: banner !< AH3 ID string
    integer :: nrOfDim_in !< corresponds to parameter nrOfDim (no. of dimensions)
    integer :: neq_in     !< corresponds to parameter neq (no. of equations)
    integer :: npr_in     !< corresponds to parameter npr (no. of processors)
    integer :: refinementFactor !< not used
    integer :: nframe           !< output counter
    real(kind=dp) :: gamma_in  !< corresponds to parameter gamma (adiab. index)
    real(kind=dp) :: time      !< output time

    ! AH3D mesh variables

    ierror=0
    
    ! Read in header
    if (rank == 0) then
       open(unit=ah3,file=filename,form="unformatted",status="old", &
            action="read")
       read(unit=ah3) banner
       read(unit=ah3) nrOfDim_in
       read(unit=ah3) neq_in
       read(unit=ah3) npr_in
       read(unit=ah3) refinementFactor
       read(unit=ah3) nframe
       read(unit=ah3) gamma_in
       read(unit=ah3) time
       
       ! Check for consistency
       if (nrOfDim_in /= nrOfDim .or. neq_in /= neq .or. npr_in /= npr .or. &
            gamma_in /= gamma ) then
          ierror=1
          write(unit=log_unit,fmt=*) &
               "Error: ah3 file inconsistent with program parameters"
       endif
       
       ! Read in mesh
       if (ierror == 0) then
          read(unit=ah3) xmesh,ymesh,zmesh
       endif

       close(unit=ah3)

    endif
    
  end subroutine restart_mesh

  !> This routine distributes n over the number of processors<br>
  !! Input:<br>
  !! nmesh    - number of mesh cells<br>
  !! numprocs - number of processors (to distribute over)<br>
  !! myid     - id of current processor<br>
  !! Output:<br>
  !! startpnt - start index for current processor<br>
  !! endpnt   - end index for current processor
  subroutine MPE_DECOMP1D (nmesh, numprocs, myid, startpnt, endpnt)

    integer,intent(in)  :: nmesh !< number of mesh cells
    integer,intent(in)  :: numprocs !< number of processors (to distribute over)
    integer,intent(in)  :: myid !< id of current processor
    integer,intent(out) :: startpnt !< start index for current processor
    integer,intent(out) :: endpnt !< end index for current processor
    integer ::  nlocal !< local number of mesh points
    integer ::  deficit !< deficit if nlocal is not a factor of nmesh
    
    nlocal  = nmesh / numprocs
    startpnt = myid * nlocal + 1
    deficit = nmesh-int(nmesh/numprocs)*numprocs !mod(nmesh,numprocs)
    startpnt = startpnt + min(myid,deficit)
    if (myid < deficit) then
       nlocal = nlocal + 1
    endif
    endpnt = startpnt + nlocal - 1
    if (endpnt > nmesh .or. myid == numprocs-1) then 
       endpnt = nmesh
    endif
  end subroutine MPE_DECOMP1D
  
  !> This routine makes the decomposition of the 3D mesh into 
  !! local processor meshes
  subroutine fnd3ddecomp ( )
    
    call MPE_DECOMP1D ( meshx, dims(1), grid_struct(1), sx, ex )
    call MPE_DECOMP1D ( meshy, dims(2), grid_struct(2), sy, ey )
    call MPE_DECOMP1D ( meshz, dims(3), grid_struct(3), sz, ez )
    
  end subroutine fnd3ddecomp
  
end module mesh
