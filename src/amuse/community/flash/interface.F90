module flash_run

! Some preprocessor directives. Now only sinkparticle_exist is used
! currently. - JW

!use particle_types
!#define TREE
#define sinkparticle_exist
!#define particle_exist
!#define no_grid

! Various interfaces and variables we modify directly in Flash. - JW

use Driver_interface, ONLY : Driver_initFlash, &
    Driver_evolveFlash, Driver_finalizeFlash, Driver_getComm, &
    Driver_getSimTime, Driver_init, Driver_getMype

use Driver_data, ONLY : dr_nbegin, dr_nend, dr_dtInit, dr_tmax, &
    dr_globalMe, dr_globalNumProcs, dr_globalComm, dr_dt, dr_dtOld, dr_dtAdvect

use Grid_interface, ONLY : Grid_getBlkData, Grid_getPointData, &
    Grid_getBlkIDFromPos, Grid_getCellCoords, Grid_getBlkIndexLimits, &
    Grid_getListOfBlocks, Grid_getBlkCenterCoords, Grid_getBlkPhysicalSize, &
    Grid_getDeltas, Grid_fillGuardCells, Grid_getBlkCornerID, Grid_getBlkPtr, &
    Grid_getBlkBoundBox, Grid_releaseBlkPtr, Grid_sortParticles, &
    Grid_getSingleCellVol, Grid_getLocalNumBlks, Grid_getMaxRefinement
    
use Grid_data, ONLY : gr_eosMode

use Gravity_interface, ONLY: Gravity_accelOneBlock, Gravity_getAccelAtPoint, &
    Gravity_getPotentialAtPoint

use Particles_interface, ONLY : Particles_sinkSyncWithParticles, &
    Particles_getGlobalNum, Particles_mapFromMesh, Particles_sinkMoveParticles, &
    Particles_sinkAccelGasOnSinks, Particles_sinkAccelSinksOnGas, &
    Particles_sinkKickGas
    
use ut_interpolationInterface, ONLY: ut_polint

! This is for Kevin Olson's tree gravity.
!def TREE

#ifdef TREE    
use tree_code_common_grav
!use amr_tree_poisson, ONLY : create_particle_tree, tree_build_grav
#endif
    
use RuntimeParameters_interface, ONLY : RuntimeParameters_set, &
    RuntimeParameters_get

#ifdef sinkparticle_exist    
    use pt_sinkInterface, ONLY : pt_sinkCreateParticle, &
    pt_sinkGatherGlobal
    
    use pt_sinkSort
#endif


#ifdef sinkparticle_exist
    use Particles_sinkdata
#endif

#ifdef particle_exist
    use Particles_data
#endif


implicit none

! Requiste header files from Flash. Note we use the mpi.h compiled
! with Flash here. This assures that we don't have an MPI mismatch.

#include "Flash.h"
#include "constants.h"
#include "Particles.h"
#define MPI
#include "Flash_mpi.h"

logical, save :: restart
character(len=16), save :: data_type, data_var

#ifdef no_grid
real :: ener_var, velx_var, vely_var, velz_var
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Runtime Parameters Settings
!!!!   These must be set before Flash
!!!!   is initialized (which is when Flash
!!!!   grabs the RT parameters info.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


FUNCTION initialize_restart()
  IMPLICIT NONE
  INTEGER :: initialize_restart
  restart = .false.
  initialize_restart=0
END FUNCTION

FUNCTION get_restart(value)
  IMPLICIT NONE
  INTEGER :: get_restart
  LOGICAL :: value
  call RuntimeParameters_get('restart',value)
!  if (restart) then
!    value=.true.
!  end if
  get_restart=0
END FUNCTION

FUNCTION set_restart(value)
  IMPLICIT NONE
  INTEGER :: set_restart
  LOGICAL :: value
  call RuntimeParameters_set('restart',value)
  restart=value
  set_restart=0
END FUNCTION
  
FUNCTION get_begin_iter_step(value)
  IMPLICIT NONE
  INTEGER :: value
  INTEGER :: get_begin_iter_step
  value = dr_nbegin
  get_begin_iter_step=0
END FUNCTION

FUNCTION set_begin_iter_step(value)
  IMPLICIT NONE
  INTEGER :: value
  INTEGER :: set_begin_iter_step
  dr_nbegin = value
  set_begin_iter_step=0
END FUNCTION

FUNCTION get_max_num_steps(value)
  IMPLICIT NONE
  INTEGER :: value
  INTEGER :: get_max_num_steps
  value = dr_nend
  get_max_num_steps=0
END FUNCTION

FUNCTION set_max_num_steps(value)
  IMPLICIT NONE
  INTEGER :: value
  INTEGER :: set_max_num_steps
  dr_nend = value
  set_max_num_steps=0
END FUNCTION

FUNCTION get_time(value)
  IMPLICIT NONE
  DOUBLE PRECISION :: value
  INTEGER :: get_time
  call Driver_getSimTime(value)
  get_time=0
END FUNCTION

FUNCTION get_end_time(value)
  IMPLICIT NONE
  DOUBLE PRECISION :: value
  INTEGER :: get_end_time
  value = dr_tmax
  get_end_time=0
END FUNCTION

FUNCTION get_timestep(value)
  IMPLICIT NONE
  DOUBLE PRECISION :: value
  INTEGER :: get_timestep
!  call RuntimeParameters_get("dtinit",value)
!!! Here it makes more sense to set it in Driver
!!! also, since a sim is already likely running.

  value = dr_dtAdvect
  get_timestep=0
END FUNCTION

FUNCTION set_timestep(value)
  IMPLICIT NONE
  DOUBLE PRECISION :: value
  INTEGER :: set_timestep
  !call RuntimeParameters_set("dtinit",value)
  dr_dt = value  ! This isn't working as intended currently. - JW
!!! Here it makes more sense to set it in Driver
!!! also, since a sim is already likely running.

  dr_dtInit = value
  set_timestep=0
END FUNCTION

FUNCTION set_end_time(value)
  IMPLICIT NONE
  DOUBLE PRECISION :: value
  INTEGER :: set_end_time
  call RuntimeParameters_set('tmax',value)
  dr_tmax = value
!  call Driver_init()
  set_end_time=0
END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Grid variable get/set operations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sets the internal energy of a block/grid.
! Note here I assumed you set the density properly first!
FUNCTION set_grid_energy_density(i, j, k, index_of_grid, nproc, enrho, n)
  IMPLICIT NONE
  INTEGER :: n, m, myProc
! cell indices, block index on local proc, proc #
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION :: en, rho, enrho(n)
  INTEGER :: set_grid_energy_density
  
call Driver_getMype(GLOBAL_COMM, myProc)

do m=1, n
  if (myProc == nproc(m)) then

        !call Grid_putBlkData(index_of_grid, CENTER, EINT_VAR, INTERIOR, [i,j,k],en)
        call Grid_getPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, [i(m),j(m),k(m)], rho)
        
        en = enrho(m) / rho
        
        call Grid_putPointData(index_of_grid(m), CENTER, EINT_VAR, INTERIOR, [i(m),j(m),k(m)], en)

  end if
end do

  set_grid_energy_density=0
END FUNCTION

! Gets the internal energy of a block/grid.
FUNCTION get_grid_energy_density(i, j, k, index_of_grid, nproc, enrho, n)
  IMPLICIT NONE
  INTEGER :: n, m, myProc, communicator, ierr
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION :: en, rho, enrho(n)
  INTEGER :: get_grid_energy_density

en=0.0
enrho=0.0

call Driver_getComm(GLOBAL_COMM, communicator)
call Driver_getMype(GLOBAL_COMM, myProc)
  
do m=1, n
  
  if (myProc == nproc(m)) then
  
    !call Grid_getBlkData(index_of_grid, CENTER, EINT_VAR, INTERIOR, [i,j,k], en)
    call Grid_getPointData(index_of_grid(m), CENTER, EINT_VAR, INTERIOR, [i(m),j(m),k(m)], en)
    call Grid_getPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, [i(m),j(m),k(m)], rho)
    
    enrho(n)=en*rho

  end if

end do

  if (myProc == 0) then
  
    call MPI_Reduce(MPI_IN_PLACE, enrho, n, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, communicator, ierr)
  else
  
    call MPI_Reduce(enrho, 0.0, n, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, communicator, ierr)
  end if


  get_grid_energy_density=0
END FUNCTION


FUNCTION get_grid_momentum_density(i, j, k, index_of_grid, nproc, & 
                                   rhovx, rhovy, rhovz, n)
  IMPLICIT NONE
  INTEGER :: n, m, myProc, communicator, ierr
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION :: rhovx(n), rhovy(n), rhovz(n), vx, vy, vz, rho
  INTEGER :: get_grid_momentum_density

!!! Note in Flash velocities are a  "per mass" variable, unlike most grid codes.

rhovx=0.0; rhovy=0.0; rhovz=0.0
vx=0.0; vy=0.0; vz=0.0; rho=0.0

call Driver_getMype(GLOBAL_COMM, myProc)
call Driver_getComm(GLOBAL_COMM, communicator)

do m=1, n 
  
  if (myProc == nproc(m)) then

      call Grid_getPointData(index_of_grid(m), CENTER, VELX_VAR, INTERIOR, [i(m),j(m),k(m)], vx)
      call Grid_getPointData(index_of_grid(m), CENTER, VELY_VAR, INTERIOR, [i(m),j(m),k(m)], vy)
      call Grid_getPointData(index_of_grid(m), CENTER, VELZ_VAR, INTERIOR, [i(m),j(m),k(m)], vz)
      call Grid_getPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, [i(m),j(m),k(m)], rho)
      
      print*, "vx =", vx, "rho =", rho
      print*, "vy =", vy, "rho =", rho
      
      rhovx(m) = rho*vx
      rhovy(m) = rho*vy
      rhovz(m) = rho*vz
  
  end if

end do

  if (myProc == 0) then
  
    call MPI_Reduce(MPI_IN_PLACE, rhovx, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, rhovy, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, rhovz, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
  else
  
    call MPI_Reduce(rhovx, 0.0, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(rhovy, 0.0, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(rhovz, 0.0, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
  end if
  


  get_grid_momentum_density=0
END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Note this function currently sets the velocity field NOT 
!!! the momentum density. Flash stores velocity not momentum.
!!! Work in progress!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Modified to assume you set the density correctly already.

FUNCTION set_grid_momentum_density(i, j, k, index_of_grid, nproc, &
                                   rhovx, rhovy, rhovz, n)
  IMPLICIT NONE
  INTEGER :: n, m, myProc
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION, dimension(n) :: rhovx, rhovy, rhovz
  DOUBLE PRECISION :: vx, vy, vz, rho
  INTEGER :: set_grid_momentum_density

do m=1, n

  call Driver_getMype(GLOBAL_COMM, myProc)
  
  if (myProc == nproc(m)) then
  
      call Grid_getPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, [i(m),j(m),k(m)], rho)

      vx = rhovx(m) / rho
      vy = rhovy(m) / rho
      vz = rhovz(m) / rho

      call Grid_putPointData(index_of_grid(m), CENTER, VELX_VAR, INTERIOR, [i(m),j(m),k(m)], vx)
      call Grid_putPointData(index_of_grid(m), CENTER, VELY_VAR, INTERIOR, [i(m),j(m),k(m)], vy)
      call Grid_putPointData(index_of_grid(m), CENTER, VELZ_VAR, INTERIOR, [i(m),j(m),k(m)], vz)
  
  end if
  
end do

  set_grid_momentum_density=0
END FUNCTION

FUNCTION get_grid_velocity(i, j, k, index_of_grid, nproc, &
                           vx, vy, vz, n)
  IMPLICIT NONE
  INTEGER :: n, m
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION, dimension(n) :: vx, vy, vz
  INTEGER :: get_grid_velocity, myProc, communicator, ierr

call Driver_getMype(GLOBAL_COMM, myProc)
call Driver_getComm(GLOBAL_COMM, communicator)

do m=1, n  
  
  if (myProc == nproc(m)) then

      call Grid_getPointData(index_of_grid(m), CENTER, VELX_VAR, &
                             INTERIOR, [i(m),j(m),k(m)], vx(m))
      call Grid_getPointData(index_of_grid(m), CENTER, VELX_VAR, &
                             INTERIOR, [i(m),j(m),k(m)], vy(m))
      call Grid_getPointData(index_of_grid(m), CENTER, VELX_VAR, &
                             INTERIOR, [i(m),j(m),k(m)], vy(m))

  end if
  
end do
  
if (myProc == 0) then
  
    call MPI_Reduce(MPI_IN_PLACE, vx, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, vy, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, vz, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
else
  
    call MPI_Reduce(vx, 0.0, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(vy, 0.0, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(vz, 0.0, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
end if

get_grid_velocity=0
END FUNCTION


FUNCTION set_grid_velocity(i, j, k, index_of_grid, nproc, vx, vy, vz, n)
  IMPLICIT NONE
  INTEGER :: n, m
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION, dimension(n) :: vx, vy, vz
  INTEGER :: set_grid_velocity, myProc

call Driver_getMype(GLOBAL_COMM, myProc)
  
do m=1, n
  
  if (myProc == nproc(m)) then

      call Grid_putPointData(index_of_grid(m), CENTER, VELX_VAR, & 
                             INTERIOR, [i(m),j(m),k(m)], vx(m))
      call Grid_putPointData(index_of_grid(m), CENTER, VELY_VAR, &
                             INTERIOR, [i(m),j(m),k(m)], vy(m))
      call Grid_putPointData(index_of_grid(m), CENTER, VELZ_VAR, &
                             INTERIOR, [i(m),j(m),k(m)], vz(m))
  
  end if

end do
set_grid_velocity=0
END FUNCTION


FUNCTION set_grid_density(i, j, k, index_of_grid, nproc, rho, n)
  IMPLICIT NONE
  INTEGER :: n, m, myProc 
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION :: rho(n)
  INTEGER :: set_grid_density
  
do m=1, n
  
  call Driver_getMype(GLOBAL_COMM, myProc)
  
  if (myProc == nproc(m)) then
  
    call Grid_putPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, [i(m),j(m),k(m)], rho(m))
  
  end if
  
end do
  
  set_grid_density=0
END FUNCTION


FUNCTION get_grid_density(i, j, k, index_of_grid, nproc, rho, n)
  IMPLICIT NONE
  INTEGER :: get_grid_density, n, m
  INTEGER :: myProc, communicator, ierr
  DOUBLE PRECISION :: rho(n)
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc

rho=0.0

call Driver_getMype(GLOBAL_COMM, myProc)
call Driver_getComm(GLOBAL_COMM, communicator)



do m=1, n

  if (myProc == nproc(m)) then

    call Grid_getPointData(index_of_grid(m),CENTER,DENS_VAR,INTERIOR,[i(m),j(m),k(m)],rho(m))

  end if

end do
  
if (myProc == 0) then

  call MPI_Reduce(MPI_IN_PLACE, rho, 1, MPI_DOUBLE_PRECISION, &
                MPI_SUM, 0, communicator, ierr)
else

  call MPI_Reduce(rho, rho, 1, MPI_DOUBLE_PRECISION, &
                MPI_SUM, 0, communicator, ierr)
end if

get_grid_density=0
END FUNCTION


FUNCTION set_grid_state(i, j, k, index_of_grid, nproc, rho, rhovx, rhovy, rhovz, rhoen, n)
  IMPLICIT NONE
  INTEGER :: n, m, myProc
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION, dimension(n) :: rho, rhovx, rhovy, rhovz, rhoen
  DOUBLE PRECISION :: vx, vy, vz, en
  INTEGER :: set_grid_state
  
do m=1, n
  
  call Driver_getMype(GLOBAL_COMM, myProc)
  
  if (myProc == nproc(m)) then
  
      vx = rhovx(m) / rho(m); vy = rhovy(m) / rho(m); vz = rhovz(m) / rho(m)
      en = rhoen(m) / rho(m)
      
      call Grid_putPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], rho(m))
      call Grid_putPointData(index_of_grid(m), CENTER, VELX_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], vx)
      call Grid_putPointData(index_of_grid(m), CENTER, VELY_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], vy)
      call Grid_putPointData(index_of_grid(m), CENTER, VELZ_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], vz)
      call Grid_putPointData(index_of_grid(m), CENTER, ENER_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], en)

  end if
  
end do
  
  set_grid_state=0
END FUNCTION

!!! TODO TODO TODO
FUNCTION get_grid_state(i, j, k, index_of_grid, nproc, &
                        rho, rhovx, rhovy, rhovz, rhoen, n)
  IMPLICIT NONE
  INTEGER :: n, m
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION, dimension(n) :: rho, rhovx, rhovy, rhovz, rhoen
  DOUBLE PRECISION :: vx, vy, vz, en
  INTEGER :: get_grid_state, myProc, communicator, ierr

  rho = 0.0
  rhovx=0.0; rhovy=0.0; rhovz=0.0
  rhoen=0.0

  call Driver_getMype(GLOBAL_COMM, myProc)
  call Driver_getComm(GLOBAL_COMM, communicator)


do m=1, n

  if (myProc == nproc(m)) then

      call Grid_getPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], rho(m))
      call Grid_getPointData(index_of_grid(m), CENTER, VELX_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], vx)
      call Grid_getPointData(index_of_grid(m), CENTER, VELX_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], vx)
      call Grid_getPointData(index_of_grid(m), CENTER, VELX_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], vx)
      call Grid_getPointData(index_of_grid(m), CENTER, ENER_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], en)
      
      rhovx(m)=vx*rho(m); rhovy(m)=rho(m)*vy; rhovz(m)=rho(m)*vz
      rhoen(m)=rho(m)*en

  end if
end do

  
  if (myProc == 0) then
  
    call MPI_Reduce(MPI_IN_PLACE, [rhovx, rhovy, rhovz, rhoen, rho], 5, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
  else
  
    call MPI_Reduce([rhovx, rhovy, rhovz, rhoen, rho], [0.0, 0.0, 0.0, 0.0, 0.0], 5, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
  end if
  
  get_grid_state=0
END FUNCTION

! This function returns the coords along a single dimension of the
! block. Limits is the number of cells in that dimension, and should
! be the same as nparts. Note that axis must be an array of length
! nparts for python to give back an array. Axis is i=1, j=2, k=3.

! Currently implemented for 1 proc only.

FUNCTION get_1blk_cell_coords(axis, blockID, limits, coords, nparts)
implicit none
integer :: get_1blk_cell_coords
integer :: axis, blockID, limits
integer :: nparts
real*8  :: coords(nparts)


call Grid_getCellCoords(axis,blockID,CENTER,.false.,coords, nparts)

get_1blk_cell_coords=0
END FUNCTION

!FUNCTION get_all_1axis_cell_coors(coords,nparts)
!integer :: get_all_1axis_cell_coors


!get_all_1axis_cell_coors=0
!END FUNCTION

!!! Here for evolve model we set tmax to the evolve time,
!!! intialize Flash, hijack the main evolution loop
!!! in Flash then set it for restart such that any evolve call
!!! after is a restart for the loop.

FUNCTION evolve_model(value)
  IMPLICIT NONE
  DOUBLE PRECISION :: value
  INTEGER :: evolve_model, num_procs, myID, ierr
!  LOGICAL, SAVE :: first_call=.true.
!  call MPI_Init(ierr)
!  call MPI_Comm_Size(MPI_COMM_WORLD, num_procs ,ierr)
!  call MPI_Comm_Rank(MPI_COMM_WORLD, myID)
  call RuntimeParameters_set('tmax',value)
  if (restart) then
    call Driver_initParallel()
    call Driver_init()
    call Driver_evolveFlash()
  else
    call Driver_initParallel()
    call Driver_init()
    call Driver_evolveFlash()
    restart=.true.
    call RuntimeParameters_set('restart',restart)
  end if
  
  evolve_model=0
END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Grid non-variable operations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION get_cell_volume(block, i, j, k, vol)
  IMPLICIT NONE
  INTEGER :: block, i, j, k, get_cell_volume
  REAL*8  :: vol
  
  call Grid_getSingleCellVol(block, INTERIOR, [i,j,k], vol)

get_cell_volume=0
END FUNCTION

! Get the number of processors which have blocks on them.
FUNCTION get_number_of_procs(n)
  implicit none
  integer :: n, get_number_of_procs
  ! Note MESH_COMM is a Flash defined constant.
  call Driver_getNumProcs(MESH_COMM, n)
  get_number_of_procs=0
END FUNCTION

FUNCTION get_all_local_num_grids(num_grid_array, nprocs)
  implicit none
  integer :: nprocs, get_all_local_num_grids
  integer :: communicator, myProc, ierr, i
  integer, dimension(nprocs) :: num_grid_array

  call Driver_getComm(GLOBAL_COMM, communicator)
  call Driver_getMype(GLOBAL_COMM, myProc)

  do i=0, nprocs-1
  
    if (myProc == i) &

      call Grid_getLocalNumBlks(num_grid_array(i+1))
      
  end do
  
  if (myProc == 0) then
  
    call MPI_REDUCE(MPI_IN_PLACE, num_grid_array, nprocs, MPI_INTEGER, & 
               MPI_SUM, 0, communicator, ierr)
  else
  
    call MPI_REDUCE(num_grid_array, num_grid_array, nprocs, MPI_INTEGER, &
                   MPI_SUM, 0, communicator, ierr)
  end if
  
  get_all_local_num_grids=0
END FUNCTION

FUNCTION get_number_of_grids(nproc, n)
  IMPLICIT NONE
  INTEGER :: n, local_n, nproc, myProc, ierr, communicator
  INTEGER :: get_number_of_grids
  INTEGER, DIMENSION(MAXBLOCKS) :: list_of_blocks
!  call MPI_Comm_Size(dr_globalComm, num_procs ,ierr)
!print*, "I'm in get grids"

!  call Driver_getNumProcs(dr_globalComm, num_procs)
!#ifdef MPI
local_n = 0
!print*, "I got the num_procs"
call Driver_getComm(GLOBAL_COMM, communicator)
call Driver_getMype(GLOBAL_COMM, myProc)

if (myProc == nproc) then

  call Grid_getLocalNumBlks(local_n)
  !print*, local_n

end if
!print*, "I called get list of blocks"

  call MPI_REDUCE(local_n, n, 1, MPI_INTEGER, MPI_SUM, 0, communicator,ierr)
!#else

!print*, "I reduced it here."

!  call Grid_getListOfBlocks(ALL_BLKS, list_of_blocks, n)
!#endif
  get_number_of_grids=0
END FUNCTION

FUNCTION set_data_type(type_in)
implicit none
integer :: set_data_type
character(len=16) :: type_in
data_type=type_in
set_data_type=0
END FUNCTION

FUNCTION set_data_var(var_in)
implicit none
integer :: set_data_var
character(len=16) :: var_in
data_var=var_in
set_data_var=0
END FUNCTION

! A function to get all of the data on the grid for a specific variable.
! Possibly we save the data_type and data_var somewhere else and then
! just get the variable here, so we don't have to pass giant arrays of them.
!FUNCTION get_data_all_local_blks(data_array, numcells)
!implicit none
!integer :: numcells, get_data_all_local_blks
!double precision, dimension(numcells) :: data_array
!integer :: myProc, nprocs, iproc, communicator, ierr
!integer, dimension(:), allocatable :: num_blocks
!integer :: blk, numblks


!call Driver_getNumProcs(nprocs)
!call Driver_getComm(MESH_COMM, communicator)
!call Driver_getMype(MESH_COMM, myProc)

!allocate(num_blocks(nprocs))


!! Get the number of blocks on each processor.
!do iproc=0, nprocs-1

!  if (myProc == iproc) then
  
!    call Grid_getLocalNumBlks(num_blocks(iproc+1))
    
!  end if
  
!end do
!! Exchange number of blocks on all processors to all processors.
!call MPI_AllReduce(MPI_IN_PLACE, num_blocks, nprocs, MPI_INTEGER, &
!                   MPI_SUM, communicator, ierr)

!if (numcells /= (nprocs*maxval(num_blocks)*NXB*NYB*NZB)) then
!    print*, "Size of data array passed doesn't equal the size of ", &
!            "nprocs*max(num_blocks)*NXB*NYB*NZB. Aborting."
!    get_data_all_local_blks = -1
!    return
!end if

!! Reshape the array, assuming it was sent with C ordering.
!data_array = reshape(data_array, [nprocs, maxval(num_blocks), NXB, NYB, NZB], order=[3,2,1])

!! Now fill the data array in parallel.
!do iproc = 0, nprocs-1

!    if (myProc == iproc) then
    
!    do blk=1, numblks

!        call Grid_getBlkData(blk, data_type, data_var, INTERIOR, &
!                             [1,1,1], data_array(iproc+1,blk,:,:,:))

!    end do
    
!    end if
    
!end do

!! Now flatten the array again. This will of course flatten in Fortran ordering.
!data_array = reshape(data_array, [numcells])

!! Reduce to the root processor.
!if (myProc == 0) then

!    call MPI_Reduce(MPI_IN_PLACE, data_array, numcells, MPI_DOUBLE_PRECISION, &
!                    MPI_SUM, communicator, ierr)
!else

!    call MPI_Reduce(data_array, 0, numcells, MPI_DOUBLE_PRECISION, &
!                    MPI_SUM, communicator, ierr)
!end if
!get_data_all_local_blks=0
!END FUNCTION

FUNCTION get_grid_range(nx, ny, nz, index_of_grid, nproc)
  IMPLICIT NONE
  INTEGER :: nx, ny, nz, nproc, myProc
  INTEGER :: index_of_grid, blkLimits(2,MDIM), blkLimitsGC(2,MDIM)
  INTEGER :: get_grid_range
  
  call Driver_getMype(GLOBAL_COMM, myProc)
  
  if (myProc == nproc) then
      
      call Grid_getBlkIndexLimits(index_of_grid, blkLimits, blkLimitsGC, CENTER)
      nx = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) +1  !Note +1 here bc I am
      ny = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) +1 !assuming indexing
      nz = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) +1 !starts at 1.
    !  imin = 1 !(blkLimitsGC(HIGH,IAXIS) - blkLimits(HIGH,IAXIS))/2 !Same here. -Josh
    !  jmin = 1 !(blkLimitsGC(HIGH,JAXIS) - blkLimits(HIGH,JAXIS))/2
    !  kmin = 1 !(blkLimitsGC(HIGH,KAXIS) - blkLimits(HIGH,KAXIS))/2

  end if
  get_grid_range=0
END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Notes on get_pos_of_index and get_index_of_pos!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Be aware that if you get the position of an index in a
!!! parent grid of Flash, then pass it back to get the index
!!! of the position, Flash is smart enough to give you the index
!!! of that location in the highest level refinement child grid
!!! instead of giving you the location in the parent grid you
!!! started with. This is a "feature" that took me a while to
!!! figure out. Working as intended! -Josh

FUNCTION get_index_of_position(x, y, z, i, j, k, index_of_grid, Proc_ID)
  IMPLICIT NONE
  INTEGER :: index_of_grid
  DOUBLE PRECISION :: x, y, z
  INTEGER :: i, j, k
  DOUBLE PRECISION, DIMENSION(MDIM) :: loc, local_pos, blockCenter, &
             blockSize, delta, stride, cornerID, cornerIDMax
  INTEGER :: get_index_of_position, blockID, blkLimits(2,MDIM), &
             blkLimitsGC(2,MDIM), Proc_ID, myProc, communicator
  INTEGER :: ii, ierr
             
  loc(1) = x; loc(2) = y; loc(3) = z
  i=0; j=0; k=0; local_pos=0.0; stride=0; cornerID=0; cornerIDMax=0
  blockID=0; delta=0.0; blockSize=0.0; blockCenter=0.0;
  blkLimits=0; blkLimitsGC=0;! Proc_ID=0; communicator=0
  index_of_grid = 0
  
  call Driver_getComm(GLOBAL_COMM, communicator)
  call Driver_getMype(GLOBAL_COMM, myProc)
  call Grid_getBlkIDFromPos(loc, blockID, Proc_ID, communicator)
  
  print*, "Comm = ", communicator
  print*, "MyProc = ", myProc
  print*, "Proc_ID = ", Proc_ID

  if (myProc == Proc_ID) then

      call Grid_getBlkCenterCoords(blockID, blockCenter)
    !  call Grid_getBlkCornerID(blockID, cornerID, stride, cornerIDMax)
      call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC, CENTER)
      call Grid_getBlkPhysicalSize(blockID, blockSize)
      call Grid_getDeltas(blockID, delta)

    !write(*,*) MDIM

    do ii=1, MDIM
      
    !  delta(ii) = blockSize(ii)/(blkLimits(HIGH,ii) - blkLimits(LOW,ii))
      local_pos(ii) = loc(ii) - blockCenter(ii) + blockSize(ii)/2.0
      
    end do
      
    !  write(*,*) loc, '\n', blockCenter, '\n', blockSize
    !  write(*,*) local_pos, '\n', delta
     
      i = ceiling(local_pos(1)/delta(1))

      if (MDIM .gt. 1) then
          j = ceiling(local_pos(2)/delta(2))
      else
          j = 0
      end if

      if (MDIM .gt. 2) then
          k = ceiling(local_pos(3)/delta(3))
      else
          k = 0
      end if

      index_of_grid = blockID
      
    !    write(*,*) i,j,k, index_of_grid

  end if
  
  if (myProc == 0) then
  
    call MPI_Reduce(MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, j, 1, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, k, 1, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, index_of_grid, 1, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
  
  else
  
    call MPI_Reduce(i, i, 1, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(j, j, 1, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(k, k, 1, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(index_of_grid, index_of_grid, 1, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    
  end if

  get_index_of_position=0
END FUNCTION

FUNCTION get_position_of_index(i, j, k, index_of_grid, Proc_ID, x, y, z)
  IMPLICIT NONE
  INTEGER :: i, j, k, index_of_grid, Proc_ID, myProc, ierr, communicator
  DOUBLE PRECISION :: x, y, z , loc(3)
  INTEGER :: get_position_of_index , indices(3)
  loc=0.0
  indices(1)=i; indices(2)=j; indices(3)=k;
  call Driver_getMype(GLOBAL_COMM, myProc)
  call Driver_getComm(GLOBAL_COMM, communicator)
  if (myProc == Proc_ID) then
    call Grid_getSingleCellCoords(indices, index_of_grid, CENTER, INTERIOR, loc)
  end if
  
  if (MyProc == 0) then
    call MPI_Reduce(MPI_IN_PLACE, loc, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
  else
    call MPI_Reduce(loc, loc, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
  end if
  
  x=loc(1); y=loc(2); z=loc(3)
  get_position_of_index=0
END FUNCTION

FUNCTION get_leaf_indices(dummy,ind,ret_cnt,num_of_blks,nparts) !
  implicit none
  integer :: get_leaf_indices, nparts, num_of_blks, ret_num_blks
  integer :: myProc, ierr, comm, i
  integer, dimension(nparts) :: ind, ret_ind, ret_cnt, dummy
  integer, dimension(dr_globalNumProcs) :: disp, rec_count
  dummy = 0
  disp = 0
  rec_count = 0
  ret_cnt = 0
  ret_ind = -1
  ind = -1
  call Grid_getListOfBlocks(LEAF, ind, num_of_blks)
  call Driver_getComm(GLOBAL_COMM, comm)
  call Driver_getMype(comm, myProc)
  !print*, "ind =", ind, myProc
  !print*, "numblks = ", num_of_blks, myProc
! Gather the array on the root process. Note that we require the
  ! user to pass the proper length of the final array.
  
  ! Make an array of the # of leaf grids from each processor.
  call MPI_Gather(num_of_blks, 1, MPI_INTEGER, &
                  rec_count, 1, MPI_INTEGER, &
                  0, comm, ierr)
  ret_cnt(:dr_globalNumProcs) = rec_count                
  ! Set the displacement for the incoming data based on how many
  ! particles are coming in from each processor. Note the displacement
  ! for the root process is zero, for rank 1 disp = num on root,
  ! for rank 2 disp = num on root + num on 1, etc etc.              
  
  do i=1, dr_globalNumProcs-1
  
    disp(i+1) = disp(i) + rec_count(i)
    
  end do
  !print*, "dr_globalNumProcs =", dr_globalNumProcs
  !print*, "rec_count = ", rec_count
  !print*, "disp = ", disp
  
  ! Now actually gather the leaf grids using the variable length array
  ! gather command in MPI.

  call MPI_Gatherv(ind, num_of_blks, MPI_INTEGER, &
                   ret_ind, rec_count, disp, MPI_INTEGER, &
                   0, comm, ierr)
  
  ind = ret_ind
  
  !print*, "ind =", ind, myProc
  
  !call MPI_Reduce(ret_num_blks,num_of_blks, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
  
  num_of_blks = sum(rec_count)
  !print*, "numblks =", num_of_blks, myProc
  
  get_leaf_indices=0
END FUNCTION

FUNCTION get_max_refinement(max_refine)
  IMPLICIT NONE
  INTEGER :: max_refine, get_max_refinement
!  INTEGER, PARAMETER :: mode=4 ! Mode 4 looks at the actual refinement level of the existing blocks.
!  INTEGER, PARAMETER :: scope=4 ! Makes comm=MPI_WORLD_COMM
  
  call Grid_getMaxRefinement(max_refine)
  
  get_max_refinement=0
END FUNCTION

FUNCTION get_potential(i, j, k, index_of_grid, potential)
  implicit none
  integer :: i, j, k, index_of_grid, get_potential
  real*8  :: potential

  call Grid_getPointData(index_of_grid, CENTER, GPOT_VAR, INTERIOR, [i, j, k], potential)
  
  get_potential=0
END FUNCTION get_potential

!!! Currently implementing!!! -JW 12-15-14
!!! Seems to be working properly now !!!  -JW 12-19-14

!FUNCTION get_potential_at_point(eps, x, y, z, potential)
!implicit none
!integer   :: get_potential_at_point
!integer   :: ProcID, communicator, blockID
!integer, parameter :: n_attrib=1
!real*8, dimension(LOW:HIGH,MDIM) :: bndbox
!integer, dimension(2,n_attrib) :: attrib
!real*8    :: eps, x, y, z
!real*8, dimension(n_attrib) :: potential
!real*8, dimension(MDIM)  :: loc, deltas
!real*8, pointer, dimension(:,:,:,:) :: solndata

!  loc = (/x, y, z/)
!  attrib(1,1) = GPOT_PART_PROP
!  attrib(2,1) = GPOT_VAR

!!print*, "Getting potential at ", loc

!  call Driver_getComm(GLOBAL_COMM, communicator)
!  call Grid_getBlkIDFromPos(loc, blockID, ProcID, communicator)
!  call Grid_getBlkPtr(blockID,solndata,CENTER)
!  call Grid_getBlkBoundBox(blockID,bndbox)
!  call Grid_getDeltas(blockID,deltas)
!  call Particles_mapFromMesh(1, 1, attrib, loc, bndbox, &
!                               deltas, solndata, potential)
!  call Grid_releaseBlkPtr(blockID,solndata)
  

!get_potential_at_point=0
!END FUNCTION

!!! This version of get_gravity matches the "spirit" of AMUSE's
!!! idea of this call. However since we want to do as few calls
!!! as possible in our production code, we use a verison that
!!! also updates the particle locations in Flash, then uses
!!! the tree to calculate this by direct sum (more acc than
!!! the finite differencing used below).

FUNCTION get_gravity_at_point(eps, x, y, z, gax, gay, gaz, nparts)
implicit none
integer   :: get_gravity_at_point, i, nparts, blkLimits(2,MDIM), &
             blkLimitsGC(2,MDIM), ii, jj, kk, iii
integer   :: ProcID, communicator, blockID, myProc, ierr, sts
integer, parameter :: n_attrib=3
real*8, dimension(LOW:HIGH,MDIM) :: bndbox
integer, dimension(2,n_attrib) :: attrib
real*8, dimension(nparts)   :: x, y, z, gax, gay, gaz
real*8    :: eps, error
real*8, dimension(nparts,MDIM) :: gravity
real*8, dimension(MDIM)  :: loc, deltas, blockSize, blockCenter, local_pos
real*8, dimension(MDIM)  :: cell_loc 
real*8, dimension(MDIM,2):: grav_cell, cell_locs
real*8, dimension(2)     :: x_cell, y_cell, z_cell
real*8, dimension(MDIM, GRID_IHI_GC, GRID_JHI_GC, GRID_KHI_GC) :: gvec
real*8, pointer, dimension(:,:,:,:) :: solndata
integer, save :: numcalls=0

!!!!! NOTE: This function REQUIRES the user to define the variable for
!!!!! graviational acceleration in each direction on the grid as
!!!!! GACX, GACY, and GACZ in the Config file. Note this is done if using
!!!!! BHTree gravity and you do bhtreeAcc=1 on setup line. See Flash 
!!!!! users guide for how to define variables for other gravity solvers.

!!!!! ALSO NOTE: The Multigrid gravity solver will give you nonsense back
!!!!! if the number of guard cells is less than four. I found out the
!!!!! hard way. Don't be like me.

!!!!! FINAL NOTE: You MUST include the unit Particles/ParticlesMapping/
!!!!! Quadratic in your setup for this function to work. If you don't
!!!!! this function will return zero gravity!

!print*, "Now we're in get_gravity_at_point."

call Grid_fillGuardCells(CENTER, ALLDIR, minLayers=NGUARD, &
                         eosMode=gr_eosMode, doEos=.false.)
                     
attrib(1,1) = ACCX_PART_PROP
attrib(2,1) = GACX_VAR
attrib(1,2) = ACCY_PART_PROP
attrib(2,2) = GACY_VAR
attrib(1,3) = ACCZ_PART_PROP
attrib(2,3) = GACZ_VAR

gravity = 0.0

 do  i=1, nparts
  
  

  loc = (/x(i), y(i), z(i)/)

!  numcalls = numcalls + 1

!  print*, "Call # ", numcalls
!  print*, "Getting gravity at ", loc

  call Driver_getMype(GLOBAL_COMM, myProc)
  call Driver_getComm(GLOBAL_COMM, communicator)
  call Grid_getBlkIDFromPos(loc, blockID, ProcID, communicator)
  
  if (myProc .eq. ProcID) then

!!!! This method uses finite differencing on the grid potential.
  
!    print*, "blockID = ", blockID

!!!! This line likely only needs to be called for the Multigrid solver.    
!    call Gravity_accelOneBlock(blockID,NGUARD,gvec)


!!    call Grid_getBlkPtr(blockID,solndata,CENTER)
    
!!    solndata(GRAX_VAR,:,:,:) = gvec(1,:,:,:)
!!    solndata(GRAY_VAR,:,:,:) = gvec(2,:,:,:)
!!    solndata(GRAZ_VAR,:,:,:) = gvec(3,:,:,:)

!!!! This method assumes we already calculated g accel when we
!!!! called Grid_solvePoisson sometime earlier. Note we might have to
!!!! now call Grid_solvePoisson during driver_init. - JW
    
    call Grid_getBlkPtr(blockID,solndata,CENTER)    
    call Grid_getBlkBoundBox(blockID,bndbox)
    call Grid_getDeltas(blockID,deltas)
    
         
!!    print*, "Cell gravity for block is: ", solndata(GACX_VAR,:,:,:) &
!!            , solndata(GACY_VAR,:,:,:), solndata(GACZ_VAR,:,:,:)
    
    call Particles_mapFromMesh(1, 3, attrib, loc, bndbox, &
                                 deltas, solndata, gravity(i,:))

    
!!    Grid_mapMeshToParticles (particles, part_props,&
!!                                    numParticles,posAttrib,&
!!                                    numAttrib, attrib,&
!!                                    mapType,gridDataStruct)
    
    
    call Grid_releaseBlkPtr(blockID,solndata)
    
!    call Grid_getBlkCenterCoords(blockID, blockCenter)
!    call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC, CENTER)
!    call Grid_getBlkPhysicalSize(blockID, blockSize)
!    call Grid_getDeltas(blockID, deltas)
    
    
    
 


!    ! Position in the block.
!    do iii=1, MDIM
!      local_pos(iii) = loc(iii) - blockCenter(iii) + blockSize(iii)/2.0
!    end do
    
!  ! Find the cell the loc is in. We use the local pos in the block including
!  ! guard cells.
   
!    ii = ceiling(local_pos(1)/deltas(1))  &
!       + (blkLimits(LOW,IAXIS)-blkLimitsGC(LOW,IAXIS))

!    if (MDIM .gt. 1) then
!        jj = ceiling(local_pos(2)/deltas(2))  &
!           + (blkLimits(LOW,JAXIS)-blkLimitsGC(LOW,JAXIS))
!    else
!        jj = 0
!    end if

!    if (MDIM .gt. 2) then
!        kk = ceiling(local_pos(3)/deltas(3))  &
!           + (blkLimits(LOW,KAXIS)-blkLimitsGC(LOW,KAXIS))
!    else
!        kk = 0
!    end if
    
!    ! Now we determine if the point is to the left or right of the
!    ! cell center of the cell it resides in. Then we get the cell centers
!    ! to the left and right of the point for linear interpolation and
!    ! the gravity to the left and right of the point in each
!    ! dimension.
    
!    if (nint(mod(local_pos(1),deltas(1))) .eq. 1) then

!        call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!                                      EXTERIOR, cell_loc)
!        cell_locs(1,1) = cell_loc(1)
!        call Grid_getSingleCellCoords([ii+1, jj, kk],blockID,CENTER, &
!                                      EXTERIOR, cell_loc)
!        cell_locs(1,2) = cell_loc(1)
         
!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!                               [ii, jj, kk],grav_cell(1,1))
!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!                               [ii+1, jj, kk],grav_cell(1,2))
!    else
    
!        call Grid_getSingleCellCoords([ii-1, jj, kk],blockID,CENTER, &
!                                      EXTERIOR, cell_loc)
!        cell_locs(1,1) = cell_loc(1)
!        call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!                                      EXTERIOR, cell_loc)
!        cell_locs(1,2) = cell_loc(1)
!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!                               [ii-1, jj, kk],grav_cell(1,1))
!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!                               [ii, jj, kk],grav_cell(1,2))
!    end if
    
!    if (nint(mod(local_pos(2),deltas(2))) .eq. 1) then

!        call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!                                      EXTERIOR, cell_loc)
!        cell_locs(2,1) = cell_loc(2)
!        call Grid_getSingleCellCoords([ii, jj+1, kk],blockID,CENTER, &
!                                      EXTERIOR, cell_loc)
!        cell_locs(2,2) = cell_loc(2)
!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!                               [ii, jj, kk],grav_cell(2,1))
!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!                               [ii, jj+1, kk],grav_cell(2,2))
!    else
    
!        call Grid_getSingleCellCoords([ii, jj-1, kk],blockID,CENTER, &
!                                      EXTERIOR, cell_loc)
!        cell_locs(2,1) = cell_loc(2)
!        call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!                                      EXTERIOR, cell_loc)
!        cell_locs(2,2) = cell_loc(2)
!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!                               [ii, jj-1, kk],grav_cell(2,1))
!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!                               [ii, jj, kk],grav_cell(2,2))
!    end if
    
!    if (nint(mod(local_pos(3),deltas(3))) .eq. 1) then

!        call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!                                      EXTERIOR, cell_loc)
!        cell_locs(3,1) = cell_loc(3)
!        call Grid_getSingleCellCoords([ii, jj, kk+1],blockID,CENTER, &
!                                  EXTERIOR, cell_loc)
!        cell_locs(3,2) = cell_loc(3)
!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!                               [ii, jj, kk],grav_cell(3,1))
!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!                               [ii, jj, kk+1],grav_cell(3,2))
!    else
    
!        call Grid_getSingleCellCoords([ii, jj, kk-1],blockID,CENTER, &
!                                      EXTERIOR, cell_loc)
!        cell_locs(3,1) = cell_loc(3)
!        call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!                                      EXTERIOR, cell_loc)
!        cell_locs(3,2) = cell_loc(3)
!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!                               [ii, jj, kk-1],grav_cell(3,1))
!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!                               [ii, jj, kk],grav_cell(3,2))
!    end if
    
!!    call Grid_getBlkPtr(blockID,solndata,CENTER)

!!    grav_cell(1,iii) = solndata(GACX_VAR,ii,jj,kk)
!!    grav_cell(2,iii) = solndata(GACY_VAR,ii,jj,kk)
!!    grav_cell(3,iii) = solndata(GACZ_VAR,ii,jj,kk)

!! Get the locations of the cells for the interpolation.
!!    call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!!                                  EXTERIOR, cell_locs(1,1))
!!    call Grid_getSingleCellCoords([ii+1, jj, kk],blockID,CENTER, &
!!                                  EXTERIOR, cell_locs(1,2))
!!    call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!!                                  EXTERIOR, cell_locs(2,1))
!!    call Grid_getSingleCellCoords([ii, jj+1, kk],blockID,CENTER, &
!!                                  EXTERIOR, cell_locs(2,2))
!!    call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!!                                  EXTERIOR, cell_locs(3,1))
!!    call Grid_getSingleCellCoords([ii, jj, kk+1],blockID,CENTER, &
!!                                  EXTERIOR, cell_locs(3,2))
                                  
!!    call Grid_releaseBlkPtr(blockID, solndata)
    
!    ! Finally we linearly interpolate the gravity between the values.
    
!    do iii=1,MDIM
    
!        call ut_polint(cell_locs(iii,:), grav_cell(iii,:), 2, &
!                       loc(iii), gravity(iii), error)
                       
!    end do
    
!    gax(i) = gravity(1); gay(i) = gravity(2); gax(i) = gravity(3)
         
!    print*, "Local cell gravity:"
!    print*, grav_cell
!    print*, "Interpolated gravity:"
!    print*, gravity
    
!    if (MyProc .ne. 0) then
    
!      call MPI_SEND(gravity, 3, MPI_DOUBLE_PRECISION, 0, 1, communicator, ierr)
      
!      print*, "Sent ", gravity
      
!    end if

  end if
    
!  if (myProc .ne. ProcID .and. myProc .eq. 0) then
  
!    call MPI_RECV(gravity, 3, MPI_DOUBLE_PRECISION, ProcID, MPI_ANY_TAG, communicator, sts, ierr)
    
!    print*, "Recieved ", gravity
  
!  end if
  
!  if (myProc .eq. 0) then
  
!    print*, "I'm zero and grav is ", gravity
    
!  end if
  
!  gax(i) = gravity(1); gay(i) = gravity(2); gaz(i) = gravity(3)
  
end do

!print*, "Out of loop."

call MPI_Reduce(gravity(:,1), gax, nparts, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
call MPI_Reduce(gravity(:,2), gay, nparts, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
call MPI_Reduce(gravity(:,3), gaz, nparts, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)

!print*, "Info sent."

!print*, "Gravity array = ", gravity, myProc
!print*, "gax =", gax, myProc

get_gravity_at_point=0
END FUNCTION


FUNCTION get_potential_at_point(eps, x, y, z, gpot, nparts)
INTEGER :: nparts, i, get_potential_at_point, MyPe, blockID, ProcID
INTEGER :: communicator, ierror
REAL*8, DIMENSION(nparts)  :: eps, x, y, z, gpot, locpot
!LOGICAL :: has_data

!has_data = .false.
locpot = 0.0
gpot = 0.0

#ifdef TREE

call Driver_getMype(GLOBAL_COMM, MyPe)
call Driver_getComm(GLOBAL_COMM, communicator)

do i=1, nparts

  call Grid_getBlkIDFromPos([x(i),y(i),z(i)], blockID, ProcID, communicator)

  if (MyPe .eq. ProcID) then

    call Gravity_getPotentialAtPoint(x(i), y(i), z(i), locpot(i))
!   has_data = .true.    
!    print*, "Local potential = ", locpot(i)

  end if

end do

!if (has_data) then

  call MPI_Reduce(locpot, gpot, nparts, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierror)

!end if

#else

do i=1, nparts

  call Gravity_getPotentialAtPoint(x(i), y(i), z(i), gpot(i))
  
end do

#endif

get_potential_at_point=0
END FUNCTION


FUNCTION get_accel_gas_on_particles(eps, x, y, z, gax, gay, gaz, nparts)
INTEGER :: nparts, i, get_accel_gas_on_particles, MyPe, blockID, ProcID
INTEGER :: communicator, ierror
REAL*8, DIMENSION(nparts)  :: eps, x, y, z, gax, gay, gaz
REAL*8, DIMENSION(nparts)  :: local_gx, local_gy, local_gz
!LOGICAL :: has_data

!has_data = .false.
local_gx=0.0; local_gy=0.0; local_gz=0.0
gax = 0.0; gay = 0.0; gaz = 0.0


! Note TREE gravity has a local tree of ALL cells on each proc
! so you only run this on the proc with the point in it.

! But multigrid only has local blocks on each proc so its run on all
! and then internally returns the summed gravity.

#ifdef TREE 
call Driver_getMype(GLOBAL_COMM, MyPe)
call Driver_getComm(GLOBAL_COMM, communicator)

do i=1, nparts

  call Grid_getBlkIDFromPos([x(i),y(i),z(i)], blockID, ProcID, communicator)

  if (MyPe .eq. ProcID) then
 
    call Gravity_getAccelAtPoint(x(i), y(i), z(i), local_gx(i), local_gy(i), local_gz(i))
!    has_data = .true.
    
  end if
    
end do



!print*, "On proc ", MyPe, "grav is", local_gx, local_gy, local_gz

!if (has_data) then
                  
!call MPI_Reduce(local_gx, gax, nparts, MPI_DOUBLE_PRECISION, &
!                  MPI_SUM, 0, communicator, ierror)
!call MPI_Reduce(local_gy, gay, nparts, MPI_DOUBLE_PRECISION, &
!                  MPI_SUM, 0, communicator, ierror)
!call MPI_Reduce(local_gz, gaz, nparts, MPI_DOUBLE_PRECISION, &
!                  MPI_SUM, 0, communicator, ierror)
!end if

#else
do i=1, nparts

    call Gravity_getAccelAtPoint(x(i), y(i), z(i), gax(i), gay(i), gaz(i))

end do
#endif

get_accel_gas_on_particles=0
END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Other Grid Maintenence. Many are not in use now.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


FUNCTION initialize_grid()
  IMPLICIT NONE
  INTEGER :: initialize_grid
  initialize_grid=0
END FUNCTION

FUNCTION initialize_code()
  IMPLICIT NONE
  INTEGER :: initialize_code
  call Driver_initParallel()
  call Driver_initFlash()
  restart = .false.
  initialize_code=0
END FUNCTION

FUNCTION cleanup_code()
  IMPLICIT NONE
  INTEGER :: cleanup_code
  call Driver_finalizeFlash()
  cleanup_code=0
END FUNCTION

FUNCTION recommit_parameters()
  IMPLICIT NONE
  INTEGER :: recommit_parameters
  recommit_parameters=0
END FUNCTION

FUNCTION commit_parameters()
  IMPLICIT NONE
  INTEGER :: commit_parameters
  commit_parameters=0
END FUNCTION

FUNCTION get_global_grid_index_limits(global_indices)
  IMPLICIT NONE
  INTEGER :: global_indices(MDIM)
  INTEGER :: get_global_grid_index_limits
  call Grid_getGlobalIndexLimits(global_indices)
  get_global_grid_index_limits=0
END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PARTICLE FUNCTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


FUNCTION get_number_of_particles(n)
  IMPLICIT NONE
  INTEGER :: n, get_number_of_particles, ierr
  n = 0
#ifdef sinkparticle_exist


call MPI_REDUCE(localnp, n, 1, MPI_INTEGER, MPI_SUM, 0, dr_globalComm,ierr)

!  n = localnp
!  n = localnpf
  
!!! Note, according to Klaus this function will sync sinks to the particle
!!! file. If needed in the future (would we run with a mix of sink and
!!! non-sink?).

!  call Particles_sinkSyncWithParticles(sink_to_part=.true.)

#endif
#ifdef particle_exist
  call Particles_getGlobalNum(n)
#endif
get_number_of_particles=0
END FUNCTION

FUNCTION get_particle_position_array(tags, x, y, z, nparts)
  implicit none
  integer :: nparts, MyPe
  double precision, dimension(nparts) :: x, y, z, tags
  integer :: get_particle_position_array, i, j, p
  integer, dimension(:), allocatable :: QSindex, id_sorted
#ifdef sinkparticle_exist

call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

! I only need one proc to do this.

if (MyPe .eq. 0) then

  ! Sort by particle tag. Note that input positions should also be
  ! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p))
  end do

  call NewQsort_IN(id_sorted, QSindex)

  ! Are we updating every particle in the simulation?
    
  if (nparts .eq. localnpf) then ! Yes, then lets do them all at once.

    do i=1, localnpf

      x(i) = particles_global(POSX_PART_PROP, QSindex(i))
      y(i) = particles_global(POSY_PART_PROP, QSindex(i))
      z(i) = particles_global(POSZ_PART_PROP, QSindex(i))
      
    end do
    
  else

  ! If not doing them all, have to do it by tag number.
    
      do j=1, nparts
      
        do i=1, localnpf
    
          if (particles_global(iptag,i) .eq. tags(j)) then
!          if (id_sorted(QSindex(i)) .eq. tags(j)) then ! Check for matching tag.

!            x(j) = particles_global(POSX_PART_PROP, QSindex(i))
!            y(j) = particles_global(POSY_PART_PROP, QSindex(i))
!            z(j) = particles_global(POSZ_PART_PROP, QSindex(i))
            x(j) = particles_global(POSX_PART_PROP, i)
            y(j) = particles_global(POSY_PART_PROP, i)
            z(j) = particles_global(POSZ_PART_PROP, i)
            
          end if
          
        end do
      
    end do
    
  end if

  deallocate(QSindex)
  deallocate(id_sorted)

end if

#endif  
get_particle_position_array=0
END FUNCTION

FUNCTION get_particle_velocity_array(tags,vx,vy,vz,nparts)
  implicit none
  integer :: nparts, MyPe
  double precision, dimension(nparts) :: vx, vy, vz, tags
  integer :: get_particle_velocity_array, i, j, p
  integer, dimension(:), allocatable :: QSindex, id_sorted
  
#ifdef sinkparticle_exist

call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

if (MyPe .eq. 0) then

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p))
  enddo

  call NewQsort_IN(id_sorted, QSindex)


  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then

  ! Yes, then lets do them all at once.

    do i=1, localnpf

      vx(i) = particles_global(VELX_PART_PROP, QSindex(i))
      vy(i) = particles_global(VELY_PART_PROP, QSindex(i))
      vz(i) = particles_global(VELZ_PART_PROP, QSindex(i))
      
    end do
    
  else

  ! If not doing them all, have to do it by tag number.
    
    do j=1, nparts
    
      do i=1, localnpf

        if (particles_global(iptag,i) .eq. tags(j)) then

!        if (id_sorted(QSindex(i)) .eq. tags(j)) then
        
!          vx(j) = particles_global(VELX_PART_PROP, QSindex(i))
!          vy(j) = particles_global(VELY_PART_PROP, QSindex(i))
!          vz(j) = particles_global(VELZ_PART_PROP, QSindex(i))
          vx(j) = particles_global(VELX_PART_PROP, i)
          vy(j) = particles_global(VELY_PART_PROP, i)
          vz(j) = particles_global(VELZ_PART_PROP, i)
          
        end if
        
      end do
      
    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)
  
end if
#endif  
get_particle_velocity_array=0
END FUNCTION

FUNCTION get_particle_mass(tags,mass,nparts)
  implicit none
  integer :: nparts, MyPe
  integer :: get_particle_mass, i, j, p
  double precision, dimension(nparts) :: mass, tags
  integer, dimension(:), allocatable :: QSindex, id_sorted
  
#ifdef sinkparticle_exist

call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

! I only need one proc to do this.

if (MyPe .eq. 0) then

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p))
  enddo

  call NewQsort_IN(id_sorted, QSindex)


  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then

  ! Yes, then lets do them all at once.

    do i=1, localnpf

      mass(i) = particles_global(ipm, QSindex(i))
      
    end do
    
  else

  ! If not doing them all, have to do it by tag number.
    
    do j=1, nparts
    
      do i=1, localnpf

        if (particles_global(iptag,i) .eq. tags(j)) then
!        if (id_sorted(QSindex(i)) .eq. int(tags(j))) then
          
          mass(j) = particles_global(ipm, i)
!          mass(j) = particles_global(ipm, QSindex(i))
          
        end if
        
      end do
      
    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)
  
end if

#endif  
get_particle_mass=0
END FUNCTION

FUNCTION set_particle_position(tags,x,y,z,nparts)
  implicit none
  integer :: nparts, MyPe
  double precision, dimension(nparts) :: x, y, z, tags
  integer :: set_particle_position, i, j, p
  integer, dimension(:), allocatable :: QSindex, id_sorted
  
#ifdef sinkparticle_exist

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p))
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?
  
if (nparts .eq. localnpf) then ! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(POSX_PART_PROP, QSindex(i)) = x(i)
    particles_global(POSY_PART_PROP, QSindex(i)) = y(i)
    particles_global(POSZ_PART_PROP, QSindex(i)) = z(i)
    
  end do
  
else

! If not doing them all, have to do it by tag number.
  
    do j=1, nparts
    
      do i=1, localnpf
  
        if (particles_global(iptag,i) .eq. tags(j)) then

!        if (id_sorted(QSindex(i)) .eq. tags(j)) then ! Check for matching tag.

!          particles_global(POSX_PART_PROP, QSindex(i)) = x(j)
!          particles_global(POSY_PART_PROP, QSindex(i)) = y(j)
!          particles_global(POSZ_PART_PROP, QSindex(i)) = z(j)
          particles_global(POSX_PART_PROP, i) = x(j)
          particles_global(POSY_PART_PROP, i) = y(j)
          particles_global(POSZ_PART_PROP, i) = z(j)
          
        end if
        
      end do
    
  end do
  
end if

deallocate(QSindex)
deallocate(id_sorted)

! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp
      
    particles_local(POSX_PART_PROP,i) = particles_global(POSX_PART_PROP, i)
    particles_local(POSY_PART_PROP,i) = particles_global(POSY_PART_PROP, i)
    particles_local(POSZ_PART_PROP,i) = particles_global(POSZ_PART_PROP, i)

end do
call pt_sinkGatherGlobal()

#endif  
set_particle_position=0
END FUNCTION

FUNCTION set_particle_velocity(tags,vx,vy,vz,nparts)
  implicit none
  integer :: nparts
  double precision, dimension(nparts) :: vx, vy, vz, tags
  integer :: set_particle_velocity, i, p, j
  integer, dimension(:), allocatable :: QSindex, id_sorted
    
#ifdef sinkparticle_exist

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p))
enddo

call NewQsort_IN(id_sorted, QSindex)


! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(VELX_PART_PROP, QSindex(i)) = vx(i)
    particles_global(VELY_PART_PROP, QSindex(i)) = vy(i)
    particles_global(VELZ_PART_PROP, QSindex(i)) = vz(i)
    
  end do
  
else

! If not doing them all, have to do it by tag number.
  
  do j=1, nparts
  
    do i=1, localnpf

       if (particles_global(iptag,i) .eq. tags(j)) then

!      if (id_sorted(QSindex(i)) .eq. tags(j)) then
      
!        particles_global(VELX_PART_PROP, QSindex(i)) = vx(j)
!        particles_global(VELY_PART_PROP, QSindex(i)) = vy(j)
!        particles_global(VELZ_PART_PROP, QSindex(i)) = vz(j)
        particles_global(VELX_PART_PROP, i) = vx(j)
        particles_global(VELY_PART_PROP, i) = vy(j)
        particles_global(VELZ_PART_PROP, i) = vz(j)
        
      end if
      
    end do
    
  end do

end if

deallocate(QSindex)
deallocate(id_sorted)

! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp
      
    particles_local(VELX_PART_PROP,i) = particles_global(VELX_PART_PROP, i)
    particles_local(VELY_PART_PROP,i) = particles_global(VELY_PART_PROP, i)
    particles_local(VELZ_PART_PROP,i) = particles_global(VELZ_PART_PROP, i)

end do
call pt_sinkGatherGlobal()

#endif  
set_particle_velocity=0
END FUNCTION

FUNCTION set_particle_mass(tags,mass, nparts)
  implicit none
  integer :: nparts
  double precision :: mass(nparts), tags(nparts)
  integer :: set_particle_mass, i, p, j, myProc, local_index, local_tag
  integer, dimension(:), allocatable :: QSindex, id_sorted

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p))
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(ipm, QSindex(i)) = mass(i)
    
  end do

else

! If not doing them all, have to do it by tag number.
  
  do j=1, nparts
  
    do i=1, localnpf
        
!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))  

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then
       
!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc
      
        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(ipm, i) = mass(j)
!        print*, particles_global(ipm, i)
        
       end if
      
    end do
    
  end do
  
end if

deallocate(QSindex)
deallocate(id_sorted)

! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp
      
    particles_local(ipm,i) = particles_global(ipm, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
call pt_sinkGatherGlobal()  
set_particle_mass=0
!stop

END FUNCTION

FUNCTION set_particle_gpot(tags,gpot,nparts)
  implicit none
  integer :: nparts
  double precision :: gpot(nparts), tags(nparts)
  integer :: set_particle_gpot, i, p, j, myProc, local_index, local_tag
  integer, dimension(:), allocatable :: QSindex, id_sorted

#ifdef sinkparticle_exist

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p))
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(GPOT_PART_PROP, QSindex(i)) = gpot(i)
    
  end do

else

! If not doing them all, have to do it by tag number.
  
  do j=1, nparts
  
    do i=1, localnpf
        
!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))  

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
      if (particles_global(iptag,i) .eq. tags(j)) then
      
!        print*, "id_sorted and global gpot", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(GPOT_PART_PROP, QSindex(i)) , myProc
!        print*, "tags and submitted gpot", tags(j), gpot(j), myProc
      
        !particles_global(GPOT_PART_PROP, QSindex(i)) = gpot(j)
        particles_global(GPOT_PART_PROP, i) = gpot(j)
        !print*, particles_global(GPOT_PART_PROP, i)
        
      end if
      
    end do
    
  end do
  
end if

deallocate(QSindex)
deallocate(id_sorted)

! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp
      
    particles_local(GPOT_PART_PROP,i) = particles_global(GPOT_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do

      
!do i=1, nparts
      
!    particles_local(GPOT_PART_PROP,n(i)) = gpot(i)

!end do
call pt_sinkGatherGlobal()
#endif  
set_particle_gpot=0
END FUNCTION

FUNCTION get_particle_gpot(n,gpot,nparts)
  implicit none
  integer :: nparts
  double precision, dimension(nparts) :: gpot
  integer :: n(nparts), get_particle_gpot, i
  
#ifdef sinkparticle_exist
      
do i=1, nparts
      
    gpot(i) = particles_local(GPOT_PART_PROP,n(i))

end do

#endif  
get_particle_gpot=0
END FUNCTION

!!! This gives the acceleration on the particles due to the gas.
!!! We will use this to add to the acceleration on the particles in
!!! PH4 during the gravity bridge.

! Note: Do we need to pass an array of the particle tags here?
! This is going to be the new get_gravity_at_point.

!!! NOTE RIGHT NOW I'M ASSUMING YOU ARE UPDATING THE PARTICLES 1 thru n
!!! IN PROPER ORDER. So x,y,z need to be in proper order.

!FUNCTION get_gravity_at_point(eps, x, y, z, gax, gay, gaz, nparts)
!!FUNCTION get_accel_gas_on_particles(eps, x, y, z, gax, gay, gaz, nparts)
!  IMPLICIT NONE
!  INTEGER :: nparts
!  DOUBLE PRECISION :: eps
!  DOUBLE PRECISION, DIMENSION(nparts) :: x, y, z, gax, gay, gaz
!!  LOGICAL :: usePart, useSinkPart
!  LOGICAL :: correct_location
!  INTEGER :: i, n(nparts) !, get_accel_gas_on_particles 
!  INTEGER :: get_gravity_at_point 
  

!#ifdef sinkparticle_exist
!!call Particles_sinkAccelGasOnSinks() ! Calculate the accel from the gas on sinks.

!call pt_sinkGatherGlobal()    

!gax=0.0; gay=0.0; gaz=0.0

!correct_location = .true.

!check_pos: do i=1, nparts

!  if (x(i) .ne. particles_local(POSX_PART_PROP, i) .or. &
!      y(i) .ne. particles_local(POSY_PART_PROP, i) .or. &
!      z(i) .ne. particles_local(POSZ_PART_PROP, i)) then
      
!      correct_location = .false.
!      exit check_pos
      
!  end if
  
!end do check_pos

!if (correct_location) then

!  call Particles_sinkAccelGasOnSinks()

!  do i=1, nparts

!!    gax(i)=0.0; gay(i)=0.0; gaz(i)=0.0
    
!    gax(i) = particles_local(ACCX_PART_PROP,i)
!    gay(i) = particles_local(ACCY_PART_PROP,i)
!    gaz(i) = particles_local(ACCZ_PART_PROP,i)

!  end do

!!  return
  
!else

!  write(*,*) "Updating particle positions."

!   n = (/ (i, i = 1, nparts) /)

!  i = set_particle_position(n, x, y, z, nparts)

!  call Particles_sinkAccelGasOnSinks()

!  do i=1, nparts

!!    gax(i)=0.0; gay(i)=0.0; gaz(i)=0.0
    
!    gax(i) = particles_local(ACCX_PART_PROP,i)
!    gay(i) = particles_local(ACCY_PART_PROP,i)
!    gaz(i) = particles_local(ACCZ_PART_PROP,i)

!  end do  

!end if

!!  else if (usePart) then
!!#elif defined particle_exist
!!    call Particles_longRangeForce()
!!    ax = particles_global(ACCX_PART_PROP,n)
!!    ay = particles_global(ACCY_PART_PROP,n)
!!    az = particles_global(ACCZ_PART_PROP,n)
!!  else
!#else
!    write(*,*) "No particles in this simulation found."
!#endif
!!  end if
!  get_gravity_at_point=0
!!  get_accel_gas_on_particles=0
!END FUNCTION

FUNCTION get_num_part_prop(n)
  implicit none
  integer :: n, get_num_part_prop
  
  n = pt_sinkParticleProps
  
  get_num_part_prop=0
END FUNCTION


FUNCTION make_sink(x, y, z, tags, nparts)
implicit none
integer :: nparts, n_created
real*8, dimension(nparts)   :: x, y, z, tags, local_tags
real*8, dimension(NPART_PROPS, nparts) :: part_copy 
integer, dimension(MAXBLOCKS, NPART_TYPES) :: particlesPerBlk
real*8   :: time
integer :: block, Proc_ID, myProc, make_sink, part_num, &
           communicator, i, ierror

! This makes the particles in parallel, but how do you know the order
! is the same as the order for AMUSE?!

n_created = 0
tags = 0.0
local_tags = 0.0

!write(*,*) nparts

do i=1, nparts

  Proc_ID = 0
  block = 0
  
  call Driver_getMype(GLOBAL_COMM, myProc)
  call Driver_getComm(GLOBAL_COMM, communicator)
  call Grid_getBlkIDFromPos([x(i),y(i),z(i)], block, Proc_ID, communicator)
  call Driver_getSimTime(time)
  
  if (myProc .eq. Proc_ID) then


    part_num = pt_sinkCreateParticle(x(i), y(i), z(i), time, block, Proc_ID)
    !n_created = n_created + 1
    local_tags(i)  = particles_local(iptag, part_num) ! should be i not n_created?
    
!    write(*,*) i, n_created
    
!  else
  
!    write(*,*) block, Proc_ID
    
  end if

end do


!write(*,*) "About to gather..."

!call MPI_BARRIER(communicator, ierror)
!call MPI_Gather(local_tags, n_created, MPI_DOUBLE_PRECISION, tags, nparts, & 
!                MPI_DOUBLE_PRECISION, communicator, 0, ierror)

call MPI_Reduce(local_tags, tags, nparts, MPI_DOUBLE_PRECISION, &
                MPI_SUM, 0, communicator, ierror)

!if (MyProc .eq. 0) then

!  call pt_sinkGatherGlobal()

!  do i=1, localnpf
    
!    print*, "Particle tag for particle ", i, " is ", particles_global(iptag,i)
!    print*, "Particle is on processor ", particles_global(ipcpu,i)
    
!  end do

!  part_copy = particles_global

!  call Grid_sortParticles(part_copy, pt_sinkParticleProps, localnpf, &
!      NPART_TYPES, sinks_maxSinks, particlesPerBlk, iptag)
      
!  do i=1, localnpf

!    print*, "Particle tag for particle ", i, " is ", part_copy(iptag,i)
!    print*, "Particle is on processor ", part_copy(ipcpu,i)
    
!  end do

!end if

call pt_sinkGatherGlobal()

!print*, "total number of particles =", localnpf

make_sink=0
END FUNCTION

FUNCTION get_particle_tags(n, tags, nparts)
integer :: get_particle_tags, nparts, i, j, p, MyPe, n(nparts)
real*8  :: tags(localnpf)
integer, dimension(:), allocatable :: QSindex, id_sorted

call pt_sinkGatherGlobal()

nparts = localnpf


call Driver_getMype(GLOBAL_COMM, MyPe)

! I only need one proc to do this.

if (MyPe .eq. 0) then

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p))
  enddo

  call NewQsort_IN(id_sorted, QSindex)

! Get all the tags. It doesn't make sense to get "some" tags,
! since tags are how we tell different particles apart.

  do i=1, localnpf

    tags(i) = particles_global(iptag, QSindex(i))
    
  end do

  deallocate(QSindex)
  deallocate(id_sorted)

end if

get_particle_tags=0

END FUNCTION

FUNCTION get_particle_proc(n, procs, nparts)
implicit none
integer :: n(nparts), get_particle_proc, nparts, i, procs(nparts), MyPe

call pt_sinkGatherGlobal()
call Driver_getMype(GLOBAL_COMM, MyPe)

if (MyPe .eq. 0) then

  do i=1, nparts

      procs(i) = particles_global(PROC_PART_PROP, n(i))
      
  end do

end if

get_particle_proc=0

END FUNCTION

FUNCTION get_particle_block(n, blocks, nparts)
implicit none
integer :: n(nparts), get_particle_block, nparts, i, blocks(nparts), MyPe

call pt_sinkGatherGlobal()
call Driver_getMype(GLOBAL_COMM, MyPe)

if (MyPe .eq. 0) then

  do i=1, nparts

      blocks(i) = particles_global(BLK_PART_PROP, n(i))
      
  end do

end if

get_particle_block=0

END FUNCTION

! Note! You have to pass an array of the number of particles for
! this function due to the interface requiring the same length
! input array as the output array. - JW

FUNCTION get_new_tags(new_tags_length, new_tags_array, nparts)
implicit none
integer   :: nparts
integer   :: new_tags_length(nparts), get_new_tags
integer   :: communicator, ierr, i
integer   :: new_tags_array(nparts), num_array(dr_globalNumProcs)
integer   :: disp(dr_globalNumProcs), MyPe, rank_minus_one
real*8    :: real_tags(nparts)

  disp = 0
  rank_minus_one = dr_globalNumProcs - 1
  call Driver_getComm(GLOBAL_COMM, communicator)
  !call Driver_getMype(GLOBAL_COMM, MyPe)
  
  ! Gather the array on the root process. Note that we require the
  ! user to pass the proper length of the final array. This can be
  ! gotten from get_number_of_new_tags.
  
  ! Make an array of the # of incoming particles from each processor.
  call MPI_Gather(number_new_sinks, 1, MPI_INTEGER, &
                  num_array, 1, MPI_INTEGER, &
                  0, communicator, ierr)
                  
  ! Set the displacement for the incoming data based on how many
  ! particles are coming in from each processor. Note the displacement
  ! for the root process is zero, for rank 1 disp = num on root,
  ! for rank 2 disp = num on root + num on 1, etc etc.              
  
  do i=1, dr_globalNumProcs-1
  
    disp(i+1) = disp(i) + num_array(i)
    
  end do
  
  ! Now actually gather the tags using the variable length array
  ! gather command in MPI.
  call MPI_Gatherv(new_tags, number_new_sinks, MPI_INTEGER, &
                  new_tags_array, num_array, disp, MPI_INTEGER, &
                  0, communicator, ierr)
                  
  !print*, "Reals are ", real_tags
                  
  !new_tags_array = floor(real_tags)
                  
  !print*, "Gathered tags are ", new_tags_array
  
  !print*, "Ierr =", ierr

!new_tags_array  = new_tags(1:new_tags_length)

get_new_tags = 0
END FUNCTION

FUNCTION get_number_of_new_tags(new_tag_num)

integer :: new_tag_num, get_number_of_new_tags
integer :: communicator, ierr

  call Driver_getComm(GLOBAL_COMM, communicator)
  
  ! Number_new_sinks is the local processor number that have
  ! been created. We combine these to give back to the interface.
  
  call MPI_Reduce(number_new_sinks, new_tag_num, 1, MPI_DOUBLE_PRECISION, &
                MPI_SUM, 0, communicator, ierr)

!new_tag_num = number_new_sinks

get_number_of_new_tags = 0
END FUNCTION

FUNCTION clear_new_tags()
implicit none
integer :: clear_new_tags

number_new_sinks = 0
new_tags = 0.0

clear_new_tags = 0
END FUNCTION

FUNCTION particles_gather()
integer :: particles_gather

call pt_sinkGatherGlobal()

particles_gather=0
END FUNCTION

!!! Give a velocity kick to one block in Flash over a time step dt.

! Here we pass accel arrays for more than one block at a time, along with
! an array with the number of blocks on each proc in ascending order (the
! same order as the accel arrays). We then loop over the accel matching
! the correct blocks with the right processor.

FUNCTION kick_block(accel_x, accel_y, accel_z, blockID, block_arr, limits, dt, nparts)
implicit none
integer :: kick_block, nparts, numBlocks, start_ind, end_ind
integer, dimension(nparts) :: blockID, block_arr
integer, dimension(dr_globalNumProcs) :: allProcs
integer :: limits, myProc, Proc_ID, communicator, ierr, ii, jj
integer, dimension(2,MDIM)  :: blkLimits, blkLimitsGC
integer :: l1,l2,l3,h1,h2,h3
real*8, dimension(nparts) :: accel_x, accel_y, accel_z
!real*8, dimension(numBlocks, limits, limits, limits) &
!                          :: acc3x, acc3y, acc3z
real*8, dimension(:,:,:,:), allocatable :: acc3x, acc3y, acc3z
real*8, dimension(limits, limits, limits) :: ovx,ovy,ovz
integer, dimension(4)     :: array_shape, order
real*8                    :: dt
real*8, dimension(MDIM)   :: blockCenter
real*8, pointer, dimension(:,:,:,:) :: solndata

call Driver_getMype(GLOBAL_COMM, myProc)
call Driver_getComm(GLOBAL_COMM, communicator)

call MPI_AllGather(myProc, 1, MPI_INTEGER, allProcs, 1, &
                MPI_INTEGER, communicator, ierr)


! Calculate the total number of blocks. This is how long the accel arrays
! have data for.
numBlocks = sum(block_arr(:dr_globalNumProcs))

!print*, "num blocks =", numBlocks

allocate(acc3x(numBlocks, limits, limits, limits))
allocate(acc3y(numBlocks, limits, limits, limits))
allocate(acc3z(numBlocks, limits, limits, limits))

! We passed the accel arrays flattened. Reform them in the proper shape
! with struct: blockID, i, j, k (cell coords).
array_shape = (/ numBlocks, limits, limits, limits /)
order = (/ 4, 3, 2, 1 /) ! Reconstruct using C ordering.
!order = (/ 1, 2, 3, 4 /)  ! Reconstruct using Fortran ordering.
! Reconstruct the arrays to match what we passed before it was flattened
! in Python.
acc3x = reshape(accel_x, array_shape, ORDER=order)
acc3y = reshape(accel_y, array_shape, ORDER=order)
acc3z = reshape(accel_z, array_shape, ORDER=order)

!print*, shape(acc3x)

!print*, acc3x(1,:,:,:)
!print*, acc3x(2,5,2,3)
!return

start_ind = 0
end_ind   = 0

! Calculate the lower and upper indices for a block on this proc.

!print*, "All procs =", allprocs

do jj=1, dr_globalNumProcs

    start_ind = 1 + end_ind
    end_ind   = sum(block_arr(:jj))

    if (myProc .eq. allProcs(jj)) then

    !print*, "start_ind = ", start_ind, myProc
    !print*, "end_in = ", end_ind, myProc

    ! Loop over all the blocks in the arrays.
        do ii=start_ind, end_ind !numBlocks

        ! Figure out which proc we are on and which one the block is on.
        ! This is clunky, but not sure how to do it better (yet).

        !call Grid_getBlkCenterCoords(blockID(ii), blockCenter)
        !call Grid_getBlkIDFromPos(blockCenter, blockID(ii), Proc_ID, communicator)

        !if (myProc .eq. Proc_ID) then

            ! Verify the limits are correct.
            call Grid_getBlkIndexLimits(blockID(ii), blkLimits, blkLimitsGC, CENTER)

            if (.NOT. (limits /= blkLimits(2,1))) then
              print*, "kick_block: Limits given don't match actual block limits. Aborting!"
              stop
            end if
            !print*, "Updating soln in block ", blockID(ii), myProc
            
            l1 = blkLimitsGC(LOW,IAXIS)+NGUARD
            h1 = blkLimitsGC(HIGH,IAXIS)-NGUARD
            l2 = blkLimitsGC(LOW,JAXIS)+NGUARD
            h2 = blkLimitsGC(HIGH,JAXIS)- NGUARD
            l3 = blkLimitsGC(LOW,KAXIS)+NGUARD
            h3 = blkLimitsGC(HIGH,KAXIS)- NGUARD
            
        !    start_ind = (ii-1)*(nparts/numBlocks)+1
        !    end_ind   = ii*(nparts/numBlocks)
            
        !    print*, "Start_ind = ", start_ind, "End_ind =", end_ind

        !    print*, acc3x(1,2,3)

            !stop            
            ! Now actually update the velocity of the gas from the kick.
            ! Note that we only kick the interior cells, not the guard cells.
            call Grid_getBlkPtr(blockID(ii),solndata,CENTER)
            
            !print*, "Random soln velx before = ", solndata(VELX_VAR,6,6,6)
            
            !print*, "acc3x =", acc3x(ii,2,2,2)*dt
            !print*, shape(acc3x(ii,:,:,:))
            !print*, shape(solndata(VELX_VAR,l1:h1,l2:h2,l3:h3))
            !print*, shape(ovx)
            !print*, "dt = ", dt
            
            !ovx = solndata(VELX_VAR,l1:h1,l2:h2,l3:h3)
            !ovy = solndata(VELY_VAR,l1:h1,l2:h2,l3:h3)
            !ovz = solndata(VELZ_VAR,l1:h1,l2:h2,l3:h3)
            
            
            solndata(VELX_VAR,l1:h1,l2:h2,l3:h3) = &
            solndata(VELX_VAR,l1:h1,l2:h2,l3:h3) + acc3x(ii,:,:,:)*dt
            solndata(VELY_VAR,l1:h1,l2:h2,l3:h3) = &
            solndata(VELY_VAR,l1:h1,l2:h2,l3:h3) + acc3y(ii,:,:,:)*dt
            solndata(VELZ_VAR,l1:h1,l2:h2,l3:h3) = &
            solndata(VELZ_VAR,l1:h1,l2:h2,l3:h3) + acc3z(ii,:,:,:)*dt

            !print*, "Random soln velx after = ", solndata(VELX_VAR,6,6,6)
            
            call Grid_releaseBlkPtr(blockID(ii),solndata)
            
        end do

    end if

end do

! This makes sure guard cells are refilled if Flash needs to do
! interpolation or averaging after we do this.

call Grid_notifySolnDataUpdate( (/VELX_VAR,VELY_VAR,VELZ_VAR/) )

deallocate(acc3x)
deallocate(acc3y)
deallocate(acc3z)

kick_block=0
END FUNCTION

FUNCTION kick_grid(dt)
implicit none
!! Give a velocity kick from gravity of sinks for time step dt.
!! ??? Then make sure the acceleration from the sinks is zeroed. ???

integer :: kick_grid, num_blks, blk_list(MAXBLOCKS), blockID
real*8  :: dt
real*8, pointer, dimension(:,:,:,:) :: solndata


call pt_sinkGatherGlobal()

call Grid_notifySolnDataUpdate( (/VELX_VAR,VELY_VAR,VELZ_VAR/) )

call Grid_getListOfBlocks(ALL_BLKS, blk_list, num_blks)

call Particles_sinkKickGas(num_blks, blk_list, dt)

!!! This is now done in Particles_sinkKickGas!!!

!do blockID=1, num_blks

!    call Grid_getBlkPtr(blockID,solndata,CENTER)

!    solndata(VELX_VAR,:,:,:) = solndata(VELX_VAR,:,:,:) + solndata(SGAX_VAR,:,:,:)*dt
!    solndata(VELY_VAR,:,:,:) = solndata(VELY_VAR,:,:,:) + solndata(SGAY_VAR,:,:,:)*dt
!    solndata(VELZ_VAR,:,:,:) = solndata(VELZ_VAR,:,:,:) + solndata(SGAZ_VAR,:,:,:)*dt
    
!    ! Zero out the accelerations so they don't get added during the
!    ! normal hydro solution.
    
!    solndata(SGAX_VAR,:,:,:) = 0.0
!    solndata(SGAY_VAR,:,:,:) = 0.0
!    solndata(SGAZ_VAR,:,:,:) = 0.0
    
!    solndata(SGXO_VAR,:,:,:) = 0.0
!    solndata(SGYO_VAR,:,:,:) = 0.0
!    solndata(SGZO_VAR,:,:,:) = 0.0
    
!    call Grid_releaseBlkPtr(blockID,solndata)
    

    
!end do

!call Grid_notifySolnDataUpdate( (/SGAX_VAR,SGAY_VAR,SGAZ_VAR,SGXO_VAR,SGYO_VAR,SGZO_VAR/) )


#ifdef TREE
if (associated(loc_t)) then
  call free_tree(loc_t)
  end if
call tree_build_grav(DENS_VAR)
#endif

kick_grid=0
END FUNCTION
#ifdef TREE
FUNCTION make_particle_tree()
implicit none
integer :: make_particle_tree

call create_particle_tree(particles_local(POSX_PART_PROP,1:localnp), &
                          particles_local(POSY_PART_PROP,1:localnp), &
                          particles_local(POSZ_PART_PROP,1:localnp), &
                          particles_local(MASS_PART_PROP,1:localnp), &
                          localnp)
                          
make_particle_tree=0
END FUNCTION
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! MORE GRID STUFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  IO stuff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION write_chpt()
implicit none
integer write_chpt

call IO_writeCheckpoint()

write_chpt=0
END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parallel initialization copied from FLASH -Josh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!****if* source/Driver/DriverMain/Driver_initParallel
!!
!! NAME
!!
!!  Driver_initParallel
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!  Initialize the parallel message-passing interface,
!!  the number of processors in a run and each processing
!!  element
!!
!!
!!  ARGUMENTS
!!    myPE : current processor
!!    numProcs : number of processors
!!  
!!
!!
!!***

!In mpif.h MPI_VERSION is an integer (and thus can't be used to conditionally
!compile code) and so we allow the user to define FLASH_MPI1 to indicate
!that they have an MPI-1 implementation.

#ifdef _OPENMP
#ifndef FLASH_MPI1
#define FLASH_MPI2_OPENMP
#endif
#endif

subroutine Driver_initParallel ()

  use Driver_data, ONLY : dr_globalMe, dr_globalNumProcs, dr_globalComm, &
       dr_mpiThreadSupport
  !$ use omp_lib
  
  implicit none             

  include "Flash_mpi.h"
  integer :: error, iprovided, errcode

#ifdef _OPENMP
#ifdef FLASH_MPI2_OPENMP
  integer, parameter :: MPI_thread_level = MPI_THREAD_SERIALIZED
#endif
#ifdef __INTEL_COMPILER
  integer(kind=kmp_size_t_kind) :: stksize
#endif
#endif
  logical :: mpiThreadSupport
  mpiThreadSupport = .false.

  !We should use MPI_Init_thread rather than MPI_Init when using multiple
  !threads so that we get a guaranteed level of thread support.

#ifdef FLASH_MPI2_OPENMP
  !We have some OpenMP parallel regions spanning MPI calls - any such
  !MPI calls are currently contained in $omp single sections and so
  !we use MPI_THREAD_SERIALIZED to give us exactly the thread support we need
  !to operate safely.  I print a warning message to the screen when your
  !MPI installation is not providing this level of thread support - it
  !is up to you whether you are happy with this risk.

  !Support Levels                     Description
  !MPI_THREAD_SINGLE     Only one thread will execute.
  !MPI_THREAD_FUNNELED   Process may be multi-threaded, but only main
  !                      thread will make MPI calls (calls are funneled to
  !                      main thread). "Default"
  !MPI_THREAD_SERIALIZED Process may be multi-threaded, any thread can
  !                      make MPI calls, but threads cannot execute MPI
  !                      calls concurrently (MPI calls are serialized).
  !MPI_THREAD_MULTIPLE   Multiple threads may call MPI, no restrictions.
   
  !The MPI standard says that "a call to MPI_INIT has the same effect as
  !a call to MPI_INIT_THREAD with a required = MPI_THREAD_SINGLE".

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These are commented out by Josh since AMUSE already loaded this! - Josh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  call MPI_Init_thread(MPI_thread_level, iprovided, error)
!  if (error /= MPI_SUCCESS) then
!     print *, "Error from MPI_Init_thread"
!     stop
!  end if
!#else
!  call MPI_Init (error)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dr_globalComm=FLASH_COMM
  call MPI_Comm_Rank (dr_globalComm, dr_globalMe, error)
  call MPI_Comm_Size (dr_globalComm, dr_globalNumProcs, error)


#ifdef _OPENMP
  if (dr_globalMe == 0) then

# ifdef FLASH_MPI2_OPENMP
     !The default thread support in Open-MPI (in the versions I have used) is
     !MPI_THREAD_SINGLE unless you configure Open-MPI with --enable-mpi-threads.

     !On Cray systems the MPI environment is limited to MPI_THREAD_SINGLE
     !by default.  This can be changed with the environmental variable
     !MPICH_MAX_THREAD_SAFETY - it has possible values of "single", "funneled",
     !"serialized" or "multiple".  To obtain MPI_THREAD_MULTIPLE thread level:
     !1) Set MPICH_MAX_THREAD_SAFETY to multiple in job submission script:
     !   export MPICH_MAX_THREAD_SAFETY="multiple"
     !2) link FLASH against a special MPI library:
     !   -lmpich_threadm.
     write(6,'(a,i3,a,i3)') " [Driver_initParallel]: "//&
          "Called MPI_Init_thread - requested level ", MPI_thread_level, &
          ", given level ", iprovided
     mpiThreadSupport = (iprovided >= MPI_thread_level);
# endif

     if (.not.mpiThreadSupport) then
        write(6,"(/ a /)") " [Driver_initParallel]: WARNING! We do not have "//&
             "a safe level of MPI thread support! (see Driver_initParallel.F90)"
        !write(6,*) "[Driver_initParalllel]: ERROR! MPI thread support too limited"
        !call MPI_Abort (dr_globalComm, errcode, error)
        !stop
     end if
  end if

  !$omp parallel
  if (dr_globalMe == 0) then
     if (omp_get_thread_num() == 0) then
        write(6,'(a,i3)') " [Driver_initParallel]: "//&
             "Number of OpenMP threads in each parallel region", &
             omp_get_num_threads()

        !Add Intel compiler specific code.  It is possible to overflow the
        !stack of the spawned OpenMP threads (e.g. WD_def 3d with block list
        !threading).  The default value for intel software stack on
        !code.uchicago.edu is 4MB (it is useful to print this information).
        !I recommend increasing this to 16MB:
        !export OMP_STACKSIZE="16M".
# ifdef __INTEL_COMPILER
        stksize = kmp_get_stacksize_s() / (1024*1024)
        write(6,'(a,i8,a)') " OpenMP thread stack size:", stksize, " MB"
# endif

        !Add Absoft compiler specific code.  The same loop iteration is
        !executed by multiple threads in parallel do loops that have 1 loop
        !iteration!  This bug happens when compiling the following test
        !problem with Absoft 64-bit Pro Fortran 11.1.4 on code.uchicago.edu.
        !
        !./setup unitTest/Multipole -auto -geometry=cartesian -3d -maxblocks=1 \
        !  +newMpole +noio threadBlockList=True -nxb=64 -nyb=64 -nzb=64
        !
        ! Set lrefine_min = lrefine_max = 1 in the flash.par.
# ifdef __ABSOFT__
        print *, ""
        print *, "WARNING!!!! Absoft compiler OpenMP bug!!!!"
        print *, "A parallel do loop with 1 loop iteration will be executed incorrectly"
        print *, ""
# endif

     end if
  end if
# ifdef DEBUG_THREADING
  write(6,'(a,i3,a,i3)') " [Driver_initParallel]: MPI rank ", dr_globalMe, &
       " has a team that includes thread ", omp_get_thread_num()
# endif
  !$omp end parallel
#endif

  dr_mpiThreadSupport = mpiThreadSupport

end subroutine Driver_initParallel






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Land of Misfit Toys. DO NOT ENTER.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!! Broken, my version with linear interpolation !!!

!FUNCTION get_potential_at_point(eps, x, y, z, potential)
!  implicit none
!  integer :: get_potential_at_point
!  INTEGER :: i, j, k, blockID, blkLimits(2,MDIM), &
!             blkLimitsGC(2,MDIM), Proc_ID, communicator
!  INTEGER :: ii, jj, kk, ind(3)
!  REAL*8, DIMENSION(MDIM) :: delta, blockCenter, &
!                             local_pos, loc, pos_in_cell, blockSize
!  real*8    :: eps, x, y, z
!  real*8    :: x1(3), x2(3), x3(3)
!  real*8    :: potential, weight(3,3,3)
!  real*8    :: gpot_cell
  
!  loc(1)=x; loc(2)=y; loc(3)=z
!  x1=0.0; x2=0.0; x3=0.0
!  weight=1.0; local_pos=0.0; pos_in_cell=0.0
!  potential = 0.0
!! Find the i,j,k index of the cell for the queried point x,y,z. 
  
!  call Driver_getComm(GLOBAL_COMM, communicator)
!  call Grid_getBlkIDFromPos(loc, blockID, Proc_ID, communicator)
!  call Grid_getBlkCenterCoords(blockID, blockCenter)
!  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC, CENTER)
!  call Grid_getBlkPhysicalSize(blockID, blockSize)
!  call Grid_getDeltas(blockID, delta)
  
!! Make sure that the guard cells in the 1st layer around the whole block
!! are filled with centered values that have been checked against the EOS
!! to be self consistent. Note eosMode is set to be the default mode
!! already being used from Grid_data.

!!  call Grid_fillGuardCells(CENTER, ALLDIR, minlayers=1, eosMode=gr_eosMode, doEos=.true.) 

!do ii=1, MDIM
!  local_pos(ii) = loc(ii) - blockCenter(ii) + blockSize(ii)/2.0 ! Position in the block.
!end do

!! Index of the cell at the local position in the block, including guard cells.
 
!  i = ceiling(local_pos(1)/delta(1)) + (blkLimits(LOW,IAXIS)-blkLimitsGC(LOW,IAXIS))

!  if (MDIM .gt. 1) then
!      j = ceiling(local_pos(2)/delta(2)) + (blkLimits(LOW,JAXIS)-blkLimitsGC(LOW,JAXIS))
!  else
!      j = 0
!  end if

!  if (MDIM .gt. 2) then
!      k = ceiling(local_pos(3)/delta(3)) + (blkLimits(LOW,KAXIS)-blkLimitsGC(LOW,KAXIS))
!  else
!      k = 0
!  end if

!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Right here we need to check for parent block and then warn. -JW
!!! NOTE: if we used stride I think it will just work... delta is always
!!! determined by the highest level of refinement at that point I think.
!!!!!!!!!!!!!!!!!!!!!!!!!

!! Where the query point x,y,z is in the local cell at i,j,k where cell size
!! is normed to 1 in each dimension and the origin is in the center.

!  pos_in_cell =  mod(local_pos,delta)/delta - 0.5
  
!!  write(*,*) i,j,k,blockID
!!  write(*,*) pos_in_cell
!!  write(*,*) local_pos
!!  write(*,*) delta
  
!! Weighting towards the inner and outer boundaries of the cell on each
!! axis for every cell that borders the cell containing the point. Note that
!! if the point is closer to one side, the cells on the opposite side give no
!! weight to the interpolation. 1 is cell "to the left" of cell with point,
!! 2 is the cell with the point, and 3 is cell "to the right" of cell with point.

!  do ii=1,MDIM

!    if (pos_in_cell(ii) > 0.0) then
!      x1(ii) = 0.0
!      x2(ii) = 1.0 - pos_in_cell(ii)
!      x3(ii) = pos_in_cell(ii)
!!      write(*,*) "Right"
    
!    else if (pos_in_cell(ii) < 0.0) then
!      x1(ii) = -pos_in_cell(ii)
!      x2(ii) = 1.0 + pos_in_cell(ii)
!      x3(ii) = 0.0
!!      write(*,*) "Left"
    
!    else
!      x1(ii) = 0.0
!      x2(ii) = 1.0
!      x3(ii) = 0.0
!!      write(*,*) "Middle"
    
!    end if
  
!  end do
  
!!  write(*,*) weight

!!do ii=1, MDIM

!  weight(1,:,:) = weight(1,:,:)*x1(1)
!  weight(2,:,:) = weight(2,:,:)*x2(1)
!  weight(3,:,:) = weight(3,:,:)*x3(1)
!  weight(:,1,:) = weight(:,1,:)*x1(2)
!  weight(:,2,:) = weight(:,2,:)*x2(2)
!  weight(:,3,:) = weight(:,3,:)*x3(2)
!  weight(:,:,1) = weight(:,:,1)*x1(3)
!  weight(:,:,2) = weight(:,:,2)*x2(3)
!  weight(:,:,3) = weight(:,:,3)*x3(3)

!!end do

!write(*,*) weight(1,2,:)
!write(*,*) weight(2,2,:)
!write(*,*) weight(3,2,:)


!  do ii=1, 3
!    do jj=1, 3
!      do kk=1, 3
!        call Grid_getPointData(blockID, CENTER, GPOT_VAR, EXTERIOR, [i-2+ii,j-2+jj,k-2+kk], gpot_cell)
!        potential = potential + weight(ii,jj,kk)*gpot_cell  
!      end do
!    end do
!  end do
  
!  get_potential_at_point=0
!END FUNCTION



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get Hydro State interpolates the data
!!! at any point. Currently broken :(
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


FUNCTION get_hydro_state_at_point(x, y, z, vx, vy, vz, &
                         rho, rhovx, rhovy, rhovz, rhoen)
  IMPLICIT NONE
  REAL*8 :: x, y, z, vx, vy, vz
  REAL*8 :: x1(3), x2(3), x3(3)
  REAL*8, DIMENSION(MDIM) :: delta, blockCenter, &
                             local_pos, loc, pos_in_cell, blockSize
  REAL*8 :: weight(3,3,3)
  REAL*8 :: rho, rhovx, rhovy, rhovz, rhoen
  REAL*8 :: rho_cell, rhovx_cell, rhovy_cell, rhovz_cell, en_cell
  INTEGER :: get_hydro_state_at_point
  INTEGER :: i, j, k, blockID, blkLimits(2,MDIM), &
             blkLimitsGC(2,MDIM), Proc_ID, communicator
  INTEGER :: ii, jj, kk
  rho=0.0; rhovx=0.0; rhovy=0.0; rhovz=0.0; rhoen = 0.0
  loc(1)=x; loc(2)=y; loc(3)=z
  x1=0.0; x2=0.0; x3=0.0
  weight=1.0; local_pos=0.0; pos_in_cell=0.0
  
  
! Find the i,j,k index of the cell for the queried point x,y,z. 
  
  call Driver_getComm(GLOBAL_COMM, communicator)
  call Grid_getBlkIDFromPos([x,y,z], blockID, Proc_ID, communicator)
  call Grid_getBlkCenterCoords(blockID, blockCenter)
  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC, CENTER)
  call Grid_getBlkPhysicalSize(blockID, blockSize)
  call Grid_getDeltas(blockID, delta)
  
! Make sure that the guard cells in the 1st layer around the whole block
! are filled with centered values that have been checked against the EOS
! to be self consistent. Note eosMode is set to be the default mode
! already being used from Grid_data.

  call Grid_fillGuardCells(CENTER, ALLDIR, minlayers=1, eosMode=gr_eosMode, doEos=.true.) 

do ii=1, MDIM
  local_pos(ii) = loc(ii) - blockCenter(ii) + blockSize(ii)/2.0 ! Position in the block.
end do
  
! Index of the cell at the local position in the block, including guard cells.
 
  i = ceiling(local_pos(1)/delta(1)) + (blkLimits(LOW,IAXIS)-blkLimitsGC(LOW,IAXIS))

  if (MDIM .gt. 1) then
      j = ceiling(local_pos(2)/delta(2)) + (blkLimits(LOW,JAXIS)-blkLimitsGC(LOW,JAXIS))
  else
      j = 0
  end if

  if (MDIM .gt. 2) then
      k = ceiling(local_pos(3)/delta(3)) + (blkLimits(LOW,KAXIS)-blkLimitsGC(LOW,KAXIS))
  else
      k = 0
  end if

! Where the query point x,y,z is in the local cell at i,j,k where cell size
! is normed to 1 in each dimension and the origin is in the center.

  pos_in_cell =  mod(local_pos,delta)/delta - 0.5
  
! Weighting towards the inner and outer boundaries of the cell on each
! axis for every cell that borders the cell containing the point. Note that
! if the point is closer to one side, the cells on the opposite side give no
! weight to the interpolation. 1 is cell "to the left" of cell with point,
! 2 is the cell with the point, and 3 is cell "to the right" of cell with point.

  do ii=1,MDIM

    if (pos_in_cell(ii) > 0.0) then
      x1(ii) = 0.0
      x2(ii) = 1.0 - pos_in_cell(ii)
      x3(ii) = pos_in_cell(ii)
    
    else if (pos_in_cell(ii) < 0.0) then
      x1(ii) = -pos_in_cell(ii)
      x2(ii) = 1.0 + pos_in_cell(ii)
      x3(ii) = 0.0
    
    else
      x1(ii) = 0.0
      x2(ii) = 1.0
      x3(ii) = 0.0
    
    end if
  
  end do

  weight(1,:,:) = weight(1,:,:)*x1(1)
  weight(2,:,:) = weight(2,:,:)*x2(1)
  weight(3,:,:) = weight(3,:,:)*x3(1)
  weight(:,1,:) = weight(:,1,:)*x1(2)
  weight(:,2,:) = weight(:,2,:)*x2(2)
  weight(:,3,:) = weight(:,3,:)*x3(2)
  weight(:,:,1) = weight(:,:,1)*x1(3)
  weight(:,:,2) = weight(:,:,2)*x2(3)
  weight(:,:,3) = weight(:,:,3)*x3(3)


!  do ii=1,3
!    do jj=1,3
!      do kk=1,3
    
!        weight(ii,jj,kk) = x1(ii)*x2(jj)*x3(kk)
    
!      end do
!    end do
!  end do

  
  do ii=1, 3
    do jj=1, 3
      do kk=1, 3
      
        call Grid_getPointData(blockID,CENTER,DENS_VAR,EXTERIOR,[i-2+ii,j-2+jj,k-2+kk],rho_cell)
        rho  = rho + weight(ii,jj,kk)*rho_cell
        call Grid_getPointData(blockID,CENTER,VELX_VAR,EXTERIOR,[i-2+ii,j-2+jj,k-2+kk],rhovx_cell)
        rhovx = rhovx + weight(ii,jj,kk)*rhovx_cell
        call Grid_getPointData(blockID,CENTER,VELY_VAR,EXTERIOR,[i-2+ii,j-2+jj,k-2+kk],rhovy_cell)
        rhovy = rhovy + weight(ii,jj,kk)*rhovy_cell
        call Grid_getPointData(blockID,CENTER,VELZ_VAR,EXTERIOR,[i-2+ii,j-2+jj,k-2+kk],rhovz_cell)
        rhovz = rhovz + weight(ii,jj,kk)*rhovz_cell
        call Grid_getPointData(blockID,CENTER,ENER_VAR,EXTERIOR,[i-2+ii,j-2+jj,k-2+kk],en_cell)
        rhoen = rhoen + weight(ii,jj,kk)*en_cell
        
      end do
    end do
  end do
  get_hydro_state_at_point = 0
END FUNCTION

!FUNCTION interpolate_at_point(x,y,z, gravity)

!  IMPLICIT NONE
!  REAL*8 :: x, y, z
!  REAL*8 :: x1(3), x2(3), x3(3)
!  REAL*8, DIMENSION(MDIM) :: delta, blockCenter, gravity, grav_cell &
!                             local_pos, loc, pos_in_cell, blockSize
!  REAL*8 :: weight(3,3,3)
!  INTEGER :: interpolate_at_point
!  INTEGER :: i, j, k, blockID, blkLimits(2,MDIM), &
!             blkLimitsGC(2,MDIM), Proc_ID, communicator
!  INTEGER :: ii, jj, kk
!  loc(1)=x; loc(2)=y; loc(3)=z
!  x1=0.0; x2=0.0; x3=0.0
!  weight=1.0; local_pos=0.0; pos_in_cell=0.0
  
  
!! Find the i,j,k index of the cell for the queried point x,y,z. 
  
!  call Driver_getComm(GLOBAL_COMM, communicator)
!  call Grid_getBlkIDFromPos([x,y,z], blockID, Proc_ID, communicator)
!  call Grid_getBlkCenterCoords(blockID, blockCenter)
!  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC, CENTER)
!  call Grid_getBlkPhysicalSize(blockID, blockSize)
!  call Grid_getDeltas(blockID, delta)
  
!! Make sure that the guard cells in the 1st layer around the whole block
!! are filled with centered values that have been checked against the EOS
!! to be self consistent. Note eosMode is set to be the default mode
!! already being used from Grid_data.

!  call Grid_fillGuardCells(CENTER, ALLDIR, minlayers=1, eosMode=gr_eosMode, doEos=.true.) 

!do ii=1, MDIM
!  local_pos(ii) = loc(ii) - blockCenter(ii) + blockSize(ii)/2.0 ! Position in the block.
!end do
  
!! Index of the cell at the local position in the block, including guard cells.
 
!  i = ceiling(local_pos(1)/delta(1)) + (blkLimits(LOW,IAXIS)-blkLimitsGC(LOW,IAXIS))

!  if (MDIM .gt. 1) then
!      j = ceiling(local_pos(2)/delta(2)) + (blkLimits(LOW,JAXIS)-blkLimitsGC(LOW,JAXIS))
!  else
!      j = 0
!  end if

!  if (MDIM .gt. 2) then
!      k = ceiling(local_pos(3)/delta(3)) + (blkLimits(LOW,KAXIS)-blkLimitsGC(LOW,KAXIS))
!  else
!      k = 0
!  end if

!! Where the query point x,y,z is in the local cell at i,j,k where cell size
!! is normed to 1 in each dimension and the origin is in the center.

!  pos_in_cell =  mod(local_pos,delta)/delta - 0.5
  
!! Weighting towards the inner and outer boundaries of the cell on each
!! axis for every cell that borders the cell containing the point. Note that
!! if the point is closer to one side, the cells on the opposite side give no
!! weight to the interpolation. 1 is cell "to the left" of cell with point,
!! 2 is the cell with the point, and 3 is cell "to the right" of cell with point.

!  do ii=1,MDIM

!    if (pos_in_cell(ii) > 0.0) then
!      x1(ii) = 0.0
!      x2(ii) = 1.0 - pos_in_cell(ii)
!      x3(ii) = pos_in_cell(ii)
    
!    else if (pos_in_cell(ii) < 0.0) then
!      x1(ii) = -pos_in_cell(ii)
!      x2(ii) = 1.0 + pos_in_cell(ii)
!      x3(ii) = 0.0
    
!    else
!      x1(ii) = 0.0
!      x2(ii) = 1.0
!      x3(ii) = 0.0
    
!    end if
  
!  end do

!  weight(1,:,:) = weight(1,:,:)*x1(1)
!  weight(2,:,:) = weight(2,:,:)*x2(1)
!  weight(3,:,:) = weight(3,:,:)*x3(1)
!  weight(:,1,:) = weight(:,1,:)*x1(2)
!  weight(:,2,:) = weight(:,2,:)*x2(2)
!  weight(:,3,:) = weight(:,3,:)*x3(2)
!  weight(:,:,1) = weight(:,:,1)*x1(3)
!  weight(:,:,2) = weight(:,:,2)*x2(3)
!  weight(:,:,3) = weight(:,:,3)*x3(3)


!!  do ii=1,3
!!    do jj=1,3
!!      do kk=1,3
    
!!        weight(ii,jj,kk) = x1(ii)*x2(jj)*x3(kk)
    
!!      end do
!!    end do
!!  end do

  
!  do ii=1, 3
!    do jj=1, 3
!      do kk=1, 3
      
!        call Grid_getPointData(blockID,CENTER,GRAX_VAR,EXTERIOR,[i-2+ii,j-2+jj,k-2+kk],grav_cell(1))
!        gravity(1) = gravity(1) + weight(ii,jj,kk)*grav_cell(1)
!        call Grid_getPointData(blockID,CENTER,GRAY_VAR,EXTERIOR,[i-2+ii,j-2+jj,k-2+kk],grav_cell(2))
!        gravity(2) = gravity(2) + weight(ii,jj,kk)*grav_cell(2)
!        call Grid_getPointData(blockID,CENTER,GRAZ_VAR,EXTERIOR,[i-2+ii,j-2+jj,k-2+kk],grav_cell(3))
!        gravity(3) = gravity(3) + weight(ii,jj,kk)*grav_cell(3)
        
!      end do
!    end do
!  end do

!END FUNCTION

end module flash_run
