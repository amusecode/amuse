module geometry

  ! Module for Capreole (f90)
  ! Author: Garrelt Mellema
  ! Date: 2007-10-02 (previous 2004-05-12, 2003-06-01)
  ! This module is also accepted by the F compiler (Dec 9, 2003)
  !
  ! This module contains the routines related to or dependent on
  ! the choice of the 3D coordinate system.
  !
  ! Version: cartesian coordinates (x,y,z).
  ! 
  ! Contents:
  ! absvel   : calculates absolute velocity
  ! timestep : calculates the CFL time step
  ! presfunc : calculates the pressure from the state variables
  !
  ! History:
  ! 2004-05-12: adapted for new approach without large (dynamic) arrays
  !             as arguments to timestep and presfunc. state is now a
  !             pointer, pointing to stnew or stold (see hydro.F90).
  ! 2007-10-02: clean up. Added "only's" to the use commands.

  use precision, only: dp
  use scaling, only: sctime, scvelo
  use sizes, only: RHO, RHVX, RHVY, RHVZ, EN, neq
  use mesh, only: sx,ex,sy,ey,sz,ez
  use grid, only: dx,dy,dz
  use atomic, only: gamma, gamma1
  use hydro, only: state,pressr,NEW,OLD,set_state_pointer
  use times, only: time

  implicit none

  private

  real(kind=dp),public :: maxdt=5.0e11_dp/SCTIME ! non CFL limit on time step
  real(kind=dp),private :: vmax_initial=10e-5_dp/SCVELO ! non CFL limit on time step

  public :: timestep,presfunc

contains

  function absvel(st0d) result(absvel_result)

    ! Function to calculate absolute velocity
    ! Version: cartesian 3D

    real(kind=dp) :: absvel_result
    real(kind=dp),dimension(neq) :: st0d

    if (st0d(RHO) == 0.0) write(*,*) "ERROR"
    absvel_result=(st0d(RHVX)*st0d(RHVX)+st0d(RHVY)*st0d(RHVY)+ &
         st0d(RHVZ)*st0d(RHVZ))/(st0d(RHO)*st0d(RHO))

  end function absvel

!--------------------------------------------------------------------

  subroutine presfunc(ist,ifi,jst,jfi,kst,kfi,newold,ierror)
    
    real(kind=dp),parameter :: precision=1.0d-13

    integer,intent(in) :: ist,ifi,jst,jfi,kst,kfi,newold
    integer,intent(out) :: ierror

    integer i,j,k

    ierror=0 ! nullify error variable

    ! Point state to appropriate array
    state => set_state_pointer(newold)

    do k=kst,kfi
       do j=jst,jfi
          do i=ist,ifi
             pressr(i,j,k)=state(i,j,k,EN)-0.5*state(i,j,k,RHO)* &
                  absvel(state(i,j,k,:))
             if (abs(pressr(i,j,k)/state(i,j,k,EN)) < precision) then
                write(30,"(A,2(1PE10.3,X),A,3(I4,X),A,E10.3)") &
                     "PRECISION ISSUE: ",pressr(i,j,k), &
                     state(i,j,k,EN)," at ",i,j,k," time = ",time
                ierror=1
             endif
             pressr(i,j,k)=gamma1*pressr(i,j,k)
             if (pressr(i,j,k) <= 0.0_dp) then
                !!write(30,"(A,2(1PE10.3,X),A,3(I4,X),A,E10.3)") &
                !!     "Presfunc reports negative pressure: ", &
                !!     pressr(i,j,k),state(i,j,k,EN)," at ",i,j,k, &
                !!     " time = ",time
                ierror=ierror+1
             endif
          enddo
       enddo
    enddo
    
  end subroutine presfunc

!------------------------------------------------------------------------------

  function timestep(cfl,newold) result(timestep_result)

    ! Function to calculate maximum time step from the CFL condition
    ! Version: spherical coordinates r, theta

    real(kind=dp) :: timestep_result

    real(kind=dp),intent(in) :: cfl
    integer,intent(in) :: newold

    real(kind=dp) :: vs ! sound speed
    real(kind=dp) :: vabs ! absolute velocities

    integer :: i,j,k ! grid loop counters
    real(kind=dp) :: vmax ! maximum signal velocity
    real(kind=dp) :: dcell ! cell size

    ! Point state to appropriate array
    state => set_state_pointer(newold)
    
    vmax=vmax_initial
    dcell=min(dy,dx,dz)
    
    do k=sz,ez
       do j=sy,ey
          do i=sx,ex
             vs=sqrt(max(0.0d0,gamma*pressr(i,j,k)/state(i,j,k,RHO)))
             vabs=sqrt(absvel(state(i,j,k,:)))
             vmax=max(vmax,(vs+vabs))
          enddo
       enddo
    enddo

    timestep_result=min(maxdt,cfl*dcell/vmax)

  end function timestep

end module geometry
