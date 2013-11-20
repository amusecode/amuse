!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters, ONLY: gravity_params
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
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp)::rx,ry,rz,rr,xc,yc,zc,eps,rcut,lr,lrmax,lrcut
  real(dp)::pi,v200,spin,c,fgas,dmax,mmin,mmax
  real(dp)::fcut,vcirc,vrot,pnfw,dnfw,mnfw,fc

  ! Add here, if you wish, some user-defined initial conditions
  eps =2.*boxlen*0.5**nlevelmax
  xc  =boxlen/2.0
  yc  =boxlen/2.0
  zc  =boxlen/2.0
  rcut=boxlen/2.0*0.75

  v200=gravity_params(1)
  c   =gravity_params(2)
  spin=gravity_params(3)
  fgas=gravity_params(4)
  dmax=1d0/eps/(1d0+eps)**2
  pi  =acos(-1d0)

  mmin=4d0*pi*dmax/3d0*eps**3
  mmax=4d0*pi*(log(1d0+rcut)-rcut/(1d0+rcut))
  vcirc=sqrt(4.*pi*(log(1d0+c)-c/(1d0+c))/c)
  fc   =c*(0.5-0.5/(1d0+c)**2-log(1d0+c)/(1d0+c))/(log(1d0+c)-c/(1d0+c))**2

  lrmax=log(10*boxlen)
  lrcut=log(rcut)
  do i=1,nn
     rx=x(i,1)-xc
     ry=x(i,2)-yc
     rz=x(i,3)-zc
     rr=sqrt(rx**2+ry**2+rz**2)
     fcut=exp(-(rr/rcut)**50)
     dnfw=1d0/rr/(1d0+rr)**2
     dnfw=MIN(dmax,dnfw)
     dnfw=MAX(dnfw*fcut,1d-6)
     mnfw=(log(1d0+rr)-rr/(1d0+rr))/(log(1d0+c)-c/(1d0+c))
     vrot=2d0*sqrt(2d0/fc)*spin*vcirc*c*mnfw/rr*fcut
     lr=MIN(log(rr),lrcut)
     pnfw=romberg(lr,lrcut,1d-4)+mmax/rr*1d-6

     q(i,1)=fgas*dnfw
     q(i,2)=-vrot/rr*ry
     q(i,3)=+vrot/rr*rx
     q(i,4)=0.0
     q(i,5)=fgas*pnfw
     if(metal)q(i,6)=z_ave*0.02
  enddo

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

contains
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function fffy(rint)
    implicit none
    !      Computes the integrand
    real(dp)::fffy
    real(dp)::rint,rrr,M,rho
    rrr=exp(rint)
    rho=1d0/rrr/(1d0+rrr)**2
    rho=MIN(rho,dmax)
    rho=MAX(rho*exp(-(rrr/rcut)**50),1d-6)
    M=4d0*pi*(log(1d0+rrr)-rrr/(1d0+rrr))
    fffy=M/rrr*rho
    return
  end function fffy
  !cccccccccccccccccccccccccccccccccccccccccccccccccccc
  function romberg(a,b,tol)
    implicit none
    real(dp)::romberg
    !
    !     Romberg returns the integral from a to b of f(x)dx using Romberg 
    !     integration. The method converges provided that f(x) is continuous 
    !     in (a,b). The function f must be double precision and must be 
    !     declared external in the calling routine.  
    !     tol indicates the desired relative accuracy in the integral.
    !
    integer::maxiter=16,maxj=5
    real(dp),dimension(100):: g
    real(dp)::a,b,tol,fourj
    real(dp)::h,error,gmax,g0,g1
    integer::nint,i,j,k,jmax

    h=0.5d0*(b-a)
    gmax=h*(fffy(a)+fffy(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if(.not.  (i>maxiter.or.(i>5.and.abs(error)<tol)))then
       !	Calculate next trapezoidal rule approximation to integral.
       
       g0=0.0d0
       do k=1,nint
          g0=g0+fffy(a+(k+k-1)*h)
       end do
       g0=0.5d0*g(1)+h*g0
       h=0.5d0*h
       nint=nint+nint
       jmax=min(i,maxj)
       fourj=1.0d0
       
       do j=1,jmax
          ! Use Richardson extrapolation.
          fourj=4.0d0*fourj
          g1=g0+(g0-g(j))/(fourj-1.0d0)
          g(j)=g0
          g0=g1
       enddo
       if (abs(g0).gt.tol) then
          error=1.0d0-gmax/g0
       else
          error=gmax
       end if
       gmax=g0
       g(jmax+1)=g0
       go to 10
    end if
    romberg=g0

    if (i>maxiter.and.abs(error)>tol) &
         &    write(*,*) 'Romberg failed to converge; integral, error=', &
         &    romberg,error

    

    return
  end function romberg

end subroutine condinit
