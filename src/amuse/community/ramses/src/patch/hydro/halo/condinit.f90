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
  real(dp)::xx,yy,zz,r,rr,rra,v,xc,yc,zc,eps,rho,rrmin,rmax,rrmax,tol,c
  real(dp)::fc,PI,IN,spin,factor,Om_b,v200

  ! Add here, if you wish, some user-defined initial conditions
  eps =gravity_params(1)
  xc  =gravity_params(2)
  yc  =gravity_params(3)
  zc  =gravity_params(4)
  v200=gravity_params(5)
  c   =gravity_params(6)
  spin=gravity_params(7)

  rmax=4*c
  tol=1d-2
  Om_b=0.1
  fc=c*(0.5-1.0/(2.0*(1.0+c)**2)-log(1.0+c)/(1.0+c)) / &
       & (log(1.0+c)-c/(1.0+c))**2
  PI=acos(-1.)
  IN=c*(c*c-3.*c-6.)/(2.*(1.0+c))+3.0*log(1.0+c)
  factor=spin*c*(log(1.+c)-c/(1.+c))*4./(sqrt(fc/2.)*PI*IN)

  do i=1,nn
     xx=x(i,1)-xc
     yy=x(i,2)-yc
     zz=x(i,3)-zc
     r  =sqrt(xx**2+yy**2)
     rr =sqrt(xx**2+yy**2+zz**2)
     rra=sqrt(xx**2+yy**2+zz**2+eps**2)
     if(rr.le.eps)then
        q(i,1)=1d0/(eps*(1.+eps)**2)
     else
        q(i,1)=1d0/(rr*(1.+rr)**2)
     endif
     q(i,1)=q(i,1)*Om_b
     v=factor*r
     q(i,ndim-1)=-v*yy/r
     q(i,ndim  )=+v*xx/r
     rrmin=log(rr)
     rrmax=log(rmax)
     q(i,ndim+2)=romberg(rrmin,rrmax,tol)
  enddo
  q(1:nn,ndim+1)=0.0

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
    real(dp)::rint,rrr,M

    rrr=exp(rint)
    ! Masse totale corrigee du facteur lie a l'auto-gravite du gaz
    !M=(log(1d0+rrr)-rrr/(1d0+rrr))
    if(rrr.le.eps)then
       rho=1d0/(eps*(1d0+eps)**2)*Om_b
       M=1d0/3d0/(eps*(1d0+eps)**2)*rrr**3
    else
       rho=1d0/(rrr*(1.+rrr)**2)*Om_b
       M=(log(1d0+rrr)-rrr/(1d0+rrr)-log(1d0+eps)+eps/(1d0+eps)) &
            & + 1d0/3d0/(1d0+eps)**2*eps**2
    endif
    M=M*(1d0+Om_b)
    fffy=c/(log(1d0+c)-c/(1d0+c))*M*rho/(rrr**2+eps**2)**(1.5)*rrr**2
    
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
