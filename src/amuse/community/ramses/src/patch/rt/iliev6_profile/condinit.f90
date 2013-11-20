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
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, r
  real(dp),parameter::n0=3.2, r0=0.0915 !cm-2 and kpc

  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

  ! Iliev7_profile patch -------------------------------------------------
  ! setup of a density profile for test 6 in Iliev's 
  ! comparison project, w. 100K temperature everywhere.
  ! the profile is:
  !                        /  n0               if   r<=r0
  !                n(r) = |
  !                        \  n0(r/r0)**-2     if   r>=r0
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  do i=1,nn
     r=sqrt(x(i,1)**2+x(i,2)**2+x(i,3)**2)! radius from origin in kpc
     if(r .le. r0) then 
        q(i,1) = n0/scale_nH
     else
        q(i,1) = n0 * (r0/r)**2 /scale_nH
     endif
     q(i,ndim+2) = 100./scale_T2 * q(i,1)
  end do
  ! Iliev7_profile patch end ---------------------------------------------


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
#if NVAR > NDIM + 2
  ! passive scalars
  do ivar=ndim+3,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit
