!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::i,j,id,iu,iv,iw,ip
  real(dp)::x0,lambday,ky,lambdaz,kz,rho1,rho2,p0,v0,v1,v2
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables

  id=1; iu=2; iv=3; iw=4; ip=ndim+2
  x0=x_center(1)
  lambday=0.25
  ky=2.*acos(-1.0d0)/lambday
  lambdaz=0.25
  kz=2.*acos(-1.0d0)/lambdaz
  rho1=d_region(1)
  rho2=d_region(2)
  v1=v_region(1)
  v2=v_region(2)
  v0=0.1
  p0=10.

  do i=1,nn
     if(x(i,1) < x0)then
        q(i,id)=rho1
        q(i,iu)=0.0
        if(ndim>1)q(i,iu)=v0*cos(ky*(x(i,2)-lambday/2.))*exp(+ky*(x(i,1)-x0))  
        if(ndim>1)q(i,iv)=v1
        if(ndim>2)q(i,iw)=0.0D0
        q(i,ip)=p0+rho1*gravity_params(1)*x(i,1)
     else
        q(i,id)=rho2
        q(i,iu)=0.0
        if(ndim>1)q(i,iu)=v0*cos(ky*(x(i,2)-lambday/2.))*exp(-ky*(x(i,1)-x0))
        if(ndim>1)q(i,iv)=v2
        if(ndim>2)q(i,iw)=0.0D0
        q(i,ip)=p0+rho1*gravity_params(1)*x0+rho2*gravity_params(1)*(x(i,1)-x0)
     endif
  end do

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
  ! passive scalars
  do ivar=ndim+3,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do

end subroutine condinit
