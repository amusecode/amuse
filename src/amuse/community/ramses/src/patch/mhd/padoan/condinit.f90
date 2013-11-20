!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft, 
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft, 
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  integer::i,j,id,iu,iv,iw,ip,ibx,iby,ibz
  real(dp)::lambda,k,rho1,p1,v1,b1,xx,yy,zz,theta,expz,v2,xp,yp,zp,expp,v3

  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

  
  id=1; iu=2; iv=3; iw=4; ip=5; ibx=6; iby=7; ibz=8
  lambda=1.
  k=2.*acos(-1.0d0)/lambda
  rho1=d_region(1)
  v1=v_region(1)
  p1=p_region(1)
  b1=b_region(1)
  v2=v_region(2)
  v3=v_region(3)
  
  theta=0.4

  do i=1,nn
     q(i,id)=rho1

     xx=x(i,1)-0.5
     yy=x(i,2)-0.5
     zz=x(i,3)-0.5
     xp=xx*cos(theta)+zz*sin(theta)
     yp=yy
     zp=zz*cos(theta)-xx*sin(theta)
     expz=exp(-3.*zp**2)
     expp=exp(-6.*(xp**2+yp**2+zp**2))
     q(i,iu)=+v1*cos(theta)*sin(k*yp)*expz                 -v2*sin(k*xp)*expp +2.0*v3
     q(i,iv)=-v1*sin(k*(xx*cos(theta)+zz*sin(theta)))*expz -v2*sin(k*yp)*expp 
     q(i,iw)=+v1*sin(theta)*sin(k*yp)*expz                 -v2*sin(k*zp)*expp -2.0*v3

     q(i,ip)=p1
     q(i,ibx)=-b1*sin(theta)
     q(i,nvar+1)=-b1*sin(theta)
     q(i,iby)=0.0
     q(i,nvar+2)=0.0
     q(i,ibz)=b1*cos(theta)
     q(i,nvar+3)=b1*cos(theta)
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
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
#if NDIM>1
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
#endif
#if NDIM>2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
#endif
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
!#if NDIM>1
!  u(1:nn,ndim+4)=q(1:nn,ndim+4)
!  u(1:nn,nvar+2)=q(1:nn,nvar+2)
!#endif
!#if NDIM>2
!  u(1:nn,ndim+5)=q(1:nn,ndim+5)
!  u(1:nn,nvar+3)=q(1:nn,nvar+3)
!#endif 
  ! passive scalars
  do ivar=9,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do

end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters  
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
!!  aa=1.0
!!  twopi=2d0*ACOS(-1d0)
!!  do i=1,ncell

!!     xx=x(i,1)
!!     yy=x(i,2)
!!     zz=x(i,3)

     ! ABC
!!     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
!!     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
!!     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

!!$     ! 1D advection test
!!$     vx=1.0_dp
!!$     vy=0.0_dp
!!$     vz=0.0_dp

!!$     ! Ponomarenko
!!$     xx=xx-boxlen/2.0
!!$     yy=yy-boxlen/2.0
!!$     rr=sqrt(xx**2+yy**2)
!!$     if(yy>0)then
!!$        tt=acos(xx/rr)
!!$     else
!!$        tt=-acos(xx/rr)+twopi
!!$     endif
!!$     if(rr<1.0)then
!!$        omega=0.609711
!!$        vz=0.792624
!!$     else
!!$        omega=0.0
!!$        vz=0.0
!!$     endif
!!$     vx=-sin(tt)*rr*omega
!!$     vy=+cos(tt)*rr*omega
     
!!     v(i,1)=vx
!!     v(i,2)=vy
!!     v(i,3)=vz

!!  end do


end subroutine velana
