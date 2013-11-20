!================================================================
!================================================================
!================================================================
!================================================================
subroutine rt_condinit(x,u,dx,nn)
  use amr_parameters
  use rt_parameters                                                    !RT
  implicit none
  integer ::nn                              ! Number of cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nrtvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim  )::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative RT variable vector. Conventions are here:
  ! U(i,1): N, U(i,2:ndim+1): Fx,Fy,Fy per group.
  ! U(:,:) are in user units.
  !================================================================

  ! Call built-in initial condition generator
  call rt_region_condinit(x,u,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

end subroutine rt_condinit
