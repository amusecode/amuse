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
  integer::ivar,i
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp)::xc,xr,xl,yl,yr,yc,zl,zr,zc,R0,V0,B0,P0,pi,f,rr,d,d_jeans,temp

  ! units
  ! rb=1
  ! G=1
  ! c_infinity=cs=1
  ! M=1

  pi=ACOS(-1.0d0)
  B0=0.0
  R0=0.25
  V0=1.0
  P0=1.0

  q(1:nn,:)=0.
  do i=1,nn

     xl=x(i,1)-0.5*dx-boxlen/2.0
     xr=x(i,1)+0.5*dx-boxlen/2.0
     xc=x(i,1)-boxlen/2.0
     yl=x(i,2)-0.5*dx-boxlen/2.0
     yr=x(i,2)+0.5*dx-boxlen/2.0
     yc=x(i,2)-boxlen/2.0
     zl=x(i,3)-0.5*dx-boxlen/2.0
     zr=x(i,3)+0.5*dx-boxlen/2.0
     zc=x(i,3)-boxlen/2.0

     rr=sqrt(xc**2+yc**2+zc**2)
!!$     if(rr.le.R0)then
!!$        ! for rr>2. dx, cloud density falls off as r^(-3/2)
!!$        q(i,1)=2.
!!$        q(i,1)=q(i,1)*MIN((rr/R0)**(-1.5),12.)
!!$     else
!!$        q(i,1)=0.02
!!$     endif

     q(i,1)=1.0
     q(i,2)=0.0
     q(i,3)=0.0
     q(i,4)=0.0
     q(i,5)=P0
     q(i,6)=B0
     q(i,nvar+1)=B0
     q(i,7)=0.0
     q(i,nvar+2)=0.0
     q(i,8)=0.0
     q(i,nvar+3)=0.0
     ! hack for heavy gas in middle
     !IF(rr .le. dx) THEN
     !   q(i,1)=10.
     !END IF
  end do

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
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
!!$  aa=1.0
!!$  twopi=2d0*ACOS(-1d0)
!!$  do i=1,ncell
!!$
!!$     xx=x(i,1)
!!$     yy=x(i,2)
!!$     zz=x(i,3)
!!$
!!$     ! ABC
!!$     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
!!$     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
!!$     vz=aa*(cos(twopi*xx)+sin(twopi*yy))
!!$
!!$     v(i,1)=vx
!!$     v(i,2)=vy
!!$     v(i,3)=vz
!!$
!!$  end do


end subroutine velana
