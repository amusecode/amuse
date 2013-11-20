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
  integer::ivar,i
  real(dp),dimension(1:nvector,1:nvar)::q   ! Primitive variables
  real(dp) :: pi,Mvirphu,rvir,cnfw,eta0,gammaKS,d_init,p_init
  real(dp) :: rx,ry,rz,xmass,ymass,zmass,rr
  ! Call built-in initial condition generator
  !call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........
  pi=3.141592
  Mvirphu = Mvir*1.99d33
  rvir = (3./4./pi*Mvirphu/(overdensity*1.89d-29*0.65**2))**0.333
  cnfw = 6./(Mvir/1.d14)**.2
  eta0 = 0.00676*(cnfw-6.5)**2+0.206*(cnfw-6.5)+2.48
  gammaKS = 1.15+0.01*(cnfw-6.5)
  d_init = 0.1*overdensity*1.89d-29*.65**2&
  &        *cnfw**2/3./(1.+cnfw)**2/(log(1.+cnfw)-cnfw/(1.+cnfw))&
  &        /(1.-3./eta0*(gammaKS-1.0d0)/gammaKS*cnfw/(log(1.+cnfw)&
  &        -cnfw/(1.+cnfw))&
  &        *(1.-log(1.+cnfw)/cnfw))**(1./(gammaKS-1.0d0))         
  p_init = 6.672d-8/3.*d_init*Mvirphu/rvir*eta0
  xmass=0.5*boxlen
  ymass=0.5*boxlen
  zmass=0.5*boxlen
  Mdot_BH=2.d72*d_init/(5./3.*p_init/d_init)**1.5

  do i=1,nn
     rx=0.0d0; ry=0.0d0; rz=0.0d0     
     rx=x(i,1)-xmass
#if NDIM>1
     ry=x(i,2)-ymass
#endif
#if NDIM>2
     rz=x(i,3)-zmass
#endif
     rr=sqrt(rx**2+ry**2+rz**2)
     U(i,1)=d_init&
     &*(1.-3./eta0*(gammaKS-1.0d0)/gammaKS*cnfw/(log(1.+cnfw)&
     &-cnfw/(1.+cnfw))&
     &*(1.-log(1.+cnfw*rr/rvir)/(cnfw*rr/rvir)))**(1./(gammaKS-1.0d0))
     U(i,2)=0.0d0 
     U(i,3)=0.0d0
     U(i,4)=0.0d0
     U(i,5)=p_init/d_init/(gamma-1.0d0)&
     &*(1.-3./eta0*(gammaKS-1.)/gammaKS*cnfw&
     &/(log(1.+cnfw)-cnfw/(1.+cnfw))&
     &*(1.-log(1.+cnfw*rr/rvir)/(cnfw*rr/rvir)))*U(i,1)
  end do


!!$  ! Convert primitive to conservative variables
!!$  ! density -> density
!!$  u(1:nn,1)=q(1:nn,1)
!!$  ! velocity -> momentum
!!$  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
!!$#if NDIM>1
!!$  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
!!$#endif
!!$#if NDIM>2
!!$  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
!!$#endif
!!$  ! kinetic energy
!!$  u(1:nn,ndim+2)=0.0d0
!!$  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
!!$#if NDIM>1
!!$  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
!!$#endif
!!$#if NDIM>2
!!$  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
!!$#endif
!!$  ! pressure -> total fluid energy
!!$  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
!!$  ! passive scalars
!!$  do ivar=ndim+3,nvar
!!$     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
!!$  end do

end subroutine condinit
