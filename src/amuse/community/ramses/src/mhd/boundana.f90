!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector
  use hydro_parameters, ONLY: nvar,boundary_var
  implicit none
  integer ::ibound                        ! Index of boundary region
  integer ::ncell                         ! Number of active cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, 
  ! U(i,6:8): Bleft, U(i,nvar+1:nvar+3): Bright
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  integer::ivar,i

  do ivar=1,nvar+3
     do i=1,ncell
        u(i,ivar)=boundary_var(ibound,ivar)
     end do
  end do
  
  ! Add here, if you wish, some user-defined boudary conditions
  ! ........

end subroutine boundana
