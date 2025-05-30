module evolution

  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2007-10-05 (previous 2004-05-11)

  ! This module contains the routines for the framework of the
  ! hydrodynamics calculation.

  ! Contents:
  ! evolve - basic integration frame work

  ! History:
  ! 2007-10-05: clean up, added only's to use statement, using
  !             log_unit

  use precision, only: dp
  use my_mpi
  use clocks, only: update_clocks
  use file_admin, only: log_unit
  use scaling, only: sctime
  use hydro, only: stold,stnew,state1,state2,state,NEW,OLD
  use times, only: time,frametime,dt,LastFrame
  use output, only: make_output
  use geometry, only: timestep
  use integrator, only: integrate

  implicit none

  real(kind=dp),parameter :: cfl=0.4d0  ! CFL number for time step (<1.0)

contains

  subroutine evolve ()
    
    ! This routine handles the basic integration frame work

    integer :: nstep,inegative,istop ! various control integers
    real(kind=dp)    :: dtlocal    ! processor local time step
    real(kind=dp)    :: nexttime   ! timer for output
    integer :: nframe              ! integers for output

#ifdef MPI
    integer :: mpi_ierror
#endif

    !--------------------------------------------------------------------------
    ! initialize output counter
    ! Note: this assumes that frametime was not changed
    ! in case of a restart
    nframe=nint(time/frametime)
    ! Report initial conditions
    if (nframe == 0) call make_output(nframe)

    ! Set time for next output
    nexttime=real(nframe+1,dp)*frametime
    nframe=nframe+1    ! update output counter

    nstep=0            ! timestep counter
    istop=0            ! initialize stop flag
    inegative=0        ! initialize negative density/energy flag
    
    ! Make stold and stnew equal to start with
    ! This assumes that the initialization routines
    ! worked with stold.
!    stnew=stold
    ! Integration loop
    do
       nstep=nstep+1 ! count the integration loops
       
       ! To avoid costly data copying, we flip pointers
       ! between state1 and state2.
       ! The integrate routine als does this. If any changes
       ! have to be made also check there.
       if (mod(nstep,2) == 0) then
          stold => state2
          stnew => state1
       else
          stold => state1
          stnew => state2
       endif
       
       ! Determine the time step on the local grid;
       ! the timestep function has to be supplied
       state => stold ! make sure state points to stold
       dtlocal=timestep(cfl,OLD)
       dt=dtlocal
#ifdef MPI
       ! communicate with all other processors to find
       ! the global minimal time step
       call MPI_ALLREDUCE(dtlocal,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
            MPI_COMM_NEW,mpi_ierror)
#endif
       ! Make sure that dt does not take us beyond nexttime
       dt=min(nexttime-time,dt)
       ! Integrate one time step
       ! istop will be non-zero in case of severe problems
       call integrate(nstep,istop)

       time=time+dt ! update time

       ! Report time step info
       write(log_unit,"(A,I6,3(X,1PE10.3))") "Time info: ",&
            nstep,time*sctime,dt*sctime,nexttime*sctime
       call flush(log_unit)

       ! Check the need for new output
       if (time >= nexttime) then
          state => stold
          call make_output(nframe)
          nframe=nframe+1
          nexttime=nexttime+frametime
       endif

       if (nframe > LastFrame .or. istop /= 0) exit ! end the integration loop

       ! Update clocks (not to loose precision in clock counters)
       call update_clocks ()

    enddo

    if (istop.ne.0) then
       ! Record conditions where error occured
       write(log_unit,*) "Stop triggered by correction routine"
       state => stold
       call make_output(nframe)
    else
       write(log_unit,*) "Maximum number of frames reached; nframe = ",nframe
    endif

  end subroutine evolve


  function evolve_step (nstep, nexttime)
    
    ! This routine handles the basic integration frame work
    integer, intent(inout) :: nstep
    integer :: evolve_step
    integer :: inegative,istop ! various control integers
    real(kind=dp)    :: dtlocal    ! processor local time step
    real(kind=dp)    :: nexttime   ! timer for output
    integer :: nframe              ! integers for output

#ifdef MPI
    integer :: mpi_ierror
#endif



    
    istop=0            ! initialize stop flag
    inegative=0        ! initialize negative density/energy flag
    
    ! Make stold and stnew equal to start with
    ! This assumes that the initialization routines
    ! worked with stold.
!    stnew=stold
    ! Integration loop
    
   nstep=nstep+1 ! count the integration loops
   
   ! To avoid costly data copying, we flip pointers
   ! between state1 and state2.
   ! The integrate routine als does this. If any changes
   ! have to be made also check there.
   if (mod(nstep,2) == 0) then
      stold => state2
      stnew => state1
   else
      stold => state1
      stnew => state2
   endif
   
   ! Determine the time step on the local grid;
   ! the timestep function has to be supplied
   state => stold ! make sure state points to stold
   dtlocal=timestep(cfl,OLD)
   dt=dtlocal
#ifdef MPI
   ! communicate with all other processors to find
   ! the global minimal time step
   call MPI_ALLREDUCE(dtlocal,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
        MPI_COMM_NEW,mpi_ierror)
#endif
   ! Make sure that dt does not take us beyond nexttime
   dt=min(nexttime-time,dt)
   ! Integrate one time step
   ! istop will be non-zero in case of severe problems
   call integrate(nstep,istop)

   time=time+dt ! update time
   ! Report time step info
   write(log_unit,"(A,I6,3(X,1PE10.3))") "Time info: ",&
        nstep,time*sctime,dt*sctime,nexttime*sctime
   call flush(log_unit)

  

   ! if (nframe > LastFrame .or. istop /= 0) exit ! end the integration loop

   ! Update clocks (not to loose precision in clock counters)
   call update_clocks ()
    
   evolve_step = istop
    
  end function evolve_step
end module evolution
