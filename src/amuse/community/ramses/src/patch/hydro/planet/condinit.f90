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
  real(dp)::x0,y0,rmax,rmin,hh,c0,rr,xx,yy,sigma,c2,pi
  real(dp)::omegaK,omega
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables

  id=1; iu=2; iv=3; iw=4; ip=ndim+2
  x0=0.5*boxlen
  y0=0.5*boxlen
  rmin=gravity_params(1)
  rmax=gravity_params(2)
  hh  =gravity_params(3)

  pi=acos(-1.0d0)

  do i=1,nn
     xx=x(i,1)-x0
     yy=x(i,2)-y0
     rr=sqrt(xx**2+yy**2)
     sigma=atan((rmax-rr)/hh)/pi + 0.5
     
     q(i,id)=MAX(sigma,smallr)

     ! Keplerian angular velocity
     omegaK=1.0/(rr**1.5+rmin**1.5)
      
     ! Sound speed profile
     c2=hh**2*(omegaK*rr)**2

     ! Rotational equilibrium angular velocity
     omega=sqrt(omegaK**2 &
          & +c2/gamma/rr**2*(2.-3.*omegaK*rr**(1.5)) + &
          & -c2/gamma/rr**2*(rr/hh)/(1.+((rmax-rr)/hh)**2)/pi/sigma )
     
     q(i,iu)=-omega*yy
     q(i,iv)=+omega*xx
     q(i,ip)=sigma*c2/gamma
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
