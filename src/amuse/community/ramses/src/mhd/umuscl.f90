! ---------------------------------------------------------------
!  MAG_UNSPLIT Unsplit second order Godunov integrator for
!              polytropic magnetized gas dynamics using 
!              MUSCL-HANCOCK scheme
!              with various Riemann solvers and slope limiters.
!              The sheme follows closely the paper by
!              Londrillo & Del Zanna ApJ 2000, 530, 508, 
!
!  inputs/outputs
!  uin         => (const)  input state
!  gravin      => (const)  input gravitational acceleration
!  iu1,iu2     => (const)  first and last index of input array,
!  ju1,ju2     => (const)  cell centered,    
!  ku1,ku2     => (const)  including buffer cells.
!  flux       <=  (modify) return fluxes in the 3 coord directions
!  if1,if2     => (const)  first and last index of output array,
!  jf1,jf2     => (const)  edge centered,
!  kf1,kf2     => (const)  for active cells only.
!  dx,dy,dz    => (const)  (dx,dy,dz)
!  dt          => (const)  time step
!  ngrid       => (const)  number of sub-grids
!  ndim        => (const)  number of dimensions
!
!  uin = (\rho, \rho u, \rho v, \rho w, Etot, A, B, C)
!  the hydro variable are cell-centered 
!  whereas the magnetic field B=(A,B,C) are face-centered.
!  Note that here we have 3 components for v and B whatever ndim.
!
!  This routine was written by Sebastien Fromang and Patrick Hennebelle
! ----------------------------------------------------------------
subroutine mag_unsplit(uin,gravin,flux,emfx,emfy,emfz,tmp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin 

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2   ,1:ndim)::tmp 

  ! Output electromotive force
  REAL(dp),DIMENSION(1:nvector,1:3,1:3,1:3)::emfx
  REAL(dp),DIMENSION(1:nvector,1:3,1:3,1:3)::emfy
  REAL(dp),DIMENSION(1:nvector,1:3,1:3,1:3)::emfz

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::qin 
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3),save::bf  

  ! Cell-centered slopes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::dq

  ! Face-centered slopes
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:ndim),save::dbf

  ! Face-averaged left and right state arrays
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qp
  
  ! Edge-averaged left-right and top-bottom state arrays
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3),save::qRT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3),save::qRB
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3),save::qLT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3),save::qLB

  ! Intermediate fluxes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::fx
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2   ),save::tx
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)       ,save::emf

  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  ! Translate to primative variables, compute sound speeds  
  call ctoprim(uin,qin,bf,gravin,dt,ngrid)

  ! Compute TVD slopes
  call uslope(bf,qin,dq,dbf,dx,dt,ngrid)

  ! Compute 3D traced-states in all three directions
#if NDIM==1
     call trace1d(qin   ,dq    ,qm,qp                ,dx      ,dt,ngrid)
#endif
#if NDIM==2
     call trace2d(qin,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy   ,dt,ngrid)
#endif
#if NDIM==3
     call trace3d(qin,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy,dz,dt,ngrid)
#endif

  ! Solve for 1D flux in X direction
  call cmpflxm(qm,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          if1  ,if2  ,jlo  ,jhi  ,klo  ,khi  , 2,3,4,6,7,8,fx,tx,ngrid)
  ! Save flux in output array
  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,1)=fx(l,i,j,k,ivar)*dt/dx
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,1)=tx(l,i,j,k,ivar)*dt/dx
        end do
     end do
  end do
  end do
  end do

  ! Solve for 1D flux in Y direction
#if NDIM>1
  call cmpflxm(qm,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jf1  ,jf2  ,klo  ,khi  , 3,2,4,7,6,8,fx,tx,ngrid)
  ! Save flux in output array
  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,2)=fx(l,i,j,k,ivar)*dt/dy
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,2)=tx(l,i,j,k,ivar)*dt/dy
        end do
     end do
  end do
  end do
  end do
#endif

  ! Solve for 1D flux in Z direction
#if NDIM==3
  call cmpflxm(qm,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jlo  ,jhi  ,kf1  ,kf2  , 4,2,3,8,6,7,fx,tx,ngrid)
  ! Save flux in output array
  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,3)=fx(l,i,j,k,ivar)*dt/dz
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,3)=tx(l,i,j,k,ivar)*dt/dz
        end do
     end do
  end do
  end do
  end do
#endif

#if NDIM>1
  ! Solve for EMF in Z direction
  CALL cmp_mag_flx(qRT,iu1+1,iu2+1,ju1+1,ju2+1,ku1  ,ku2  , &
       &           qRB,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &           qLT,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &           qLB,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &               if1  ,if2  ,jf1  ,jf2  ,klo  ,khi  , 2,3,4,6,7,8,emf,ngrid)
 ! Save vector in output array
  do k=klo,khi
  do j=jf1,jf2
  do i=if1,if2
     do l=1,ngrid
        emfz(l,i,j,k)=emf(l,i,j,k)*dt/dx
     end do
  end do
  end do
  end do
#if NDIM==2
  do k=2,2
  do j=jf1,jf2
  do i=if1,if2
     do l=1,ngrid
        emfz(l,i,j,k)=emf(l,i,j,k-1)*dt/dx
     end do
  end do
  end do
  end do
#endif
#endif

#if NDIM>2
  ! Solve for EMF in Y direction
  CALL cmp_mag_flx(qRT,iu1+1,iu2+1,ju1,ju2,ku1+1,ku2+1, &
       &           qLT,iu1  ,iu2  ,ju1,ju2,ku1+1,ku2+1, &
       &           qRB,iu1+1,iu2+1,ju1,ju2,ku1  ,ku2  , &
       &           qLB,iu1  ,iu2  ,ju1,ju2,ku1  ,ku2  , &
       &               if1  ,if2  ,jlo,jhi,kf1  ,kf2  , 4,2,3,8,6,7,emf,ngrid)
  ! Save vector in output array
  do k=kf1,kf2
  do j=jlo,jhi
  do i=if1,if2
     do l=1,ngrid
        emfy(l,i,j,k)=emf(l,i,j,k)*dt/dx
     end do
  end do
  end do
  end do
  ! Solve for EMF in X direction
  CALL cmp_mag_flx(qRT,iu1,iu2,ju1+1,ju2+1,ku1+1,ku2+1, &
       &           qRB,iu1,iu2,ju1+1,ju2+1,ku1  ,ku2  , &
       &           qLT,iu1,iu2,ju1  ,ju2  ,ku1+1,ku2+1, &
       &           qLB,iu1,iu2,ju1  ,ju2  ,ku1  ,ku2  , &
       &               ilo,ihi,jf1  ,jf2  ,kf1  ,kf2  , 3,4,2,7,8,6,emf,ngrid)
  ! Save vector in output array
  do k=kf1,kf2
  do j=jf1,jf2
  do i=ilo,ihi
     do l=1,ngrid
        emfx(l,i,j,k)=emf(l,i,j,k)*dt/dx
     end do
  end do
  end do
  end do
#endif

end subroutine mag_unsplit
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM==1
SUBROUTINE  trace1d(q,dq,qm,qp,dx,dt,ngrid)
  USE amr_parameters
  USE hydro_parameters
  USE const
  IMPLICIT NONE

  INTEGER ::ngrid
  REAL(dp)::dx,dt

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  ! declare local variables
  INTEGER ::i, j, k, l, n
  INTEGER ::ilo,ihi,jlo,jhi,klo,khi
  INTEGER ::ir, iu, iv, iw, ip, iA, iB, iC
  REAL(dp)::dtdx
  REAL(dp)::r, u, v, w, p, A, B, C
  REAL(dp)::drx, dux, dvx, dwx, dpx, dAx, dBx, dCx
  REAL(dp)::sr0, su0=0, sv0=0, sw0=0, sp0, sA0, sB0, sC0
  
  dtdx = dt/dx

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir = 1; iu = 2; iv = 3 ; iw = 4 ; ip = 5; iA = 6; iB = 7; iC = 8

  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi
           DO l = 1, ngrid

              ! Cell centered values
              r = q(l,i,j,k,ir)
              u = q(l,i,j,k,iu)
              v = q(l,i,j,k,iv)
              w = q(l,i,j,k,iw)
              p = q(l,i,j,k,ip)
              A = q(l,i,j,k,iA)
              B = q(l,i,j,k,iB)
              C = q(l,i,j,k,iC)

              ! TVD slopes in X direction
              drx = half*dq(l,i,j,k,ir,1)
              dux = half*dq(l,i,j,k,iu,1)
              dvx = half*dq(l,i,j,k,iv,1)
              dwx = half*dq(l,i,j,k,iw,1)
              dpx = half*dq(l,i,j,k,ip,1)
              dBx = half*dq(l,i,j,k,iB,1)
              dCx = half*dq(l,i,j,k,iC,1)

              ! Source terms (including transverse derivatives)
              sr0 = -u*drx-r*dux
              if(scheme.ne.'induction')then
              su0 = -u*dux-(dpx+B*dBx+C*dCx)/r  
              sv0 = -u*dvx+(A*dBx)/r
              sw0 = -u*dwx+(A*dCx)/r
              endif
              sp0 = -u*dpx-gamma*p*dux
              sB0 = -u*dBx+A*dvx-B*dux
              sC0 = -u*dCx+A*dwx-C*dux

              ! Cell-centered predicted states
              r = r + sr0*dtdx
              u = u + su0*dtdx
              v = v + sv0*dtdx
              w = w + sw0*dtdx
              p = p + sp0*dtdx
              B = B + sB0*dtdx
              C = C + sC0*dtdx

              ! Right state at left interface
              qp(l,i,j,k,ir,1) = r - drx
              qp(l,i,j,k,iu,1) = u - dux
              qp(l,i,j,k,iv,1) = v - dvx
              qp(l,i,j,k,iw,1) = w - dwx
              qp(l,i,j,k,ip,1) = p - dpx
              qp(l,i,j,k,iA,1) = A
              qp(l,i,j,k,iB,1) = B - dBx
              qp(l,i,j,k,iC,1) = C - dCx
              qp(l,i,j,k,ir,1) = MAX(smallr, qp(l,i,j,k,ir,1))

              ! Left state at right interface
              qm(l,i,j,k,ir,1) = r + drx
              qm(l,i,j,k,iu,1) = u + dux
              qm(l,i,j,k,iv,1) = v + dvx
              qm(l,i,j,k,iw,1) = w + dwx
              qm(l,i,j,k,ip,1) = p + dpx
              qm(l,i,j,k,iA,1) = A
              qm(l,i,j,k,iB,1) = B + dBx
              qm(l,i,j,k,iC,1) = C + dCx
              qm(l,i,j,k,ir,1) = MAX(smallr, qm(l,i,j,k,ir,1))
           END DO
        END DO
     END DO
  END DO
  
  ! passive scalars
#if NVAR > 8
  DO n = 9, nvar
     DO k = klo, khi
        DO j = jlo, jhi
           DO i = ilo, ihi
              DO l = 1, ngrid
                 a   = q(l,i,j,k,n )           ! Cell centered values
                 u   = q(l,i,j,k,iu)  
                 dax = half * dq(l,i,j,k,n,1)  ! TVD slope
                 sa0 = -u*dax                  ! Source terms
                 a   = a + sa0*dtdx            ! Predicted state
                 qp(l,i,j,k,n,1) = a - dax     ! Right state
                 qm(l,i,j,k,n,1) = a + dax     ! Left state
              END DO
           END DO
        END DO
     END DO
  END DO
#endif

END SUBROUTINE trace1d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM==2
SUBROUTINE trace2d(q,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy,dt,ngrid)
  USE amr_parameters
  USE hydro_parameters
  USE const
  IMPLICIT NONE

  INTEGER ::ngrid
  REAL(dp)::dx, dy, dt

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  REAL(dp),DIMENSION(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:ndim)::dbf

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRB
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLB

  ! Declare local variables
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::Ez
  INTEGER ::i, j, k, l, n
  INTEGER ::ilo,ihi,jlo,jhi,klo,khi
  INTEGER ::ir, iu, iv, iw, ip, iA, iB, iC 
  REAL(dp)::dtdx, dtdy, smallp
  REAL(dp)::r, u, v, w, p, A, B, C
  REAL(dp)::ELL, ELR, ERL, ERR
  REAL(dp)::drx, dux, dvx, dwx, dpx, dAx, dBx, dCx
  REAL(dp)::dry, duy, dvy, dwy, dpy, dAy, dBy, dCy
  REAL(dp)::sr0, su0=0, sv0=0, sw0=0, sp0, sA0, sB0, sC0
  REAL(dp)::AL, AR, BL, BR
  REAL(dp)::dALy, dARy, dBLx, dBRx
  REAL(DP)::sAL0, sAR0, sBL0, sBR0
  
  dtdx = dt/dx
  dtdy = dt/dy
  smallp = smallr*smallc**2/gamma

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5; ia=6; ib=7; ic=8 

  DO k = klo, ku2
     DO j = jlo, ju2
        DO i = ilo, iu2
           DO l = 1, ngrid
              u = 0.25*(q(l,i-1,j-1,k,iu)+q(l,i-1,j,k,iu)+q(l,i,j-1,k,iu)+q(l,i,j,k,iu))
              v = 0.25*(q(l,i-1,j-1,k,iv)+q(l,i-1,j,k,iv)+q(l,i,j-1,k,iv)+q(l,i,j,k,iv))
              A = 0.5*(bf(l,i,j-1,k,1)+bf(l,i,j,k,1))
              B = 0.5*(bf(l,i-1,j,k,2)+bf(l,i,j,k,2))
              Ez(l,i,j,k)=u*B-v*A
           END DO
        END DO
     END DO
  END DO

  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi
           DO l = 1, ngrid

              ! Cell centered values
              r =    q(l,i,j,k,ir)
              u =    q(l,i,j,k,iu)
              v =    q(l,i,j,k,iv)
              w =    q(l,i,j,k,iw)
              p =    q(l,i,j,k,ip)
              A =    q(l,i,j,k,iA)
              B =    q(l,i,j,k,iB)
              C =    q(l,i,j,k,iC)

              ! Face centered variables
              AL =  bf(l,i  ,j  ,k,1)
              AR =  bf(l,i+1,j  ,k,1)
              BL =  bf(l,i  ,j  ,k,2)
              BR =  bf(l,i  ,j+1,k,2)

              ! Cell centered TVD slopes in X direction
              drx = half * dq(l,i,j,k,ir,1)
              dux = half * dq(l,i,j,k,iu,1)
              dvx = half * dq(l,i,j,k,iv,1)
              dwx = half * dq(l,i,j,k,iw,1)
              dpx = half * dq(l,i,j,k,ip,1)
              dBx = half * dq(l,i,j,k,iB,1)
              dCx = half * dq(l,i,j,k,iC,1)

              ! Cell centered TVD slopes in Y direction
              dry = half * dq(l,i,j,k,ir,2)
              duy = half * dq(l,i,j,k,iu,2)
              dvy = half * dq(l,i,j,k,iv,2)
              dwy = half * dq(l,i,j,k,iw,2)
              dpy = half * dq(l,i,j,k,ip,2)
              dAy = half * dq(l,i,j,k,iA,2)
              dCy = half * dq(l,i,j,k,iC,2)

              ! Face centered TVD slopes in transverse direction
              dALy = half * dbf(l,i  ,j  ,k,1,1)
              dARy = half * dbf(l,i+1,j  ,k,1,1)
              dBLx = half * dbf(l,i  ,j  ,k,2,1)
              dBRx = half * dbf(l,i  ,j+1,k,2,1)

              ! Edge centered electric field Ez = uB-vA
              ELL = Ez(l,i  ,j  ,k)
              ELR = Ez(l,i  ,j+1,k)
              ERL = Ez(l,i+1,j  ,k)
              ERR = Ez(l,i+1,j+1,k)

              ! Face-centered predicted states
              sAL0 = +(ELR-ELL)*dtdy*half
              sAR0 = +(ERR-ERL)*dtdy*half
              sBL0 = -(ERL-ELL)*dtdx*half
              sBR0 = -(ERR-ELR)*dtdx*half

              AL = AL + sAL0
              AR = AR + sAR0
              BL = BL + sBL0
              BR = BR + sBR0
              
              ! Source terms (including transverse derivatives)
              sr0 = (-u*drx-dux*r)*dtdx + (-v*dry-dvy*r)*dtdy
              if(scheme.ne.'induction')then
              su0 = (-u*dux-(dpx+B*dBx+C*dCx)/r)*dtdx + (-v*duy+B*dAy/r)*dtdy 
              sv0 = (-u*dvx+A*dBx/r)*dtdx + (-v*dvy-(dpy+A*dAy+C*dCy)/r)*dtdy
              sw0 = (-u*dwx+A*dCx/r)*dtdx + (-v*dwy+B*dCy/r)*dtdy
              endif
              sp0 = (-u*dpx-dux*gamma*p)*dtdx + (-v*dpy-dvy*gamma*p)*dtdy
              sC0 = (-u*dCx-C*dux+A*dwx)*dtdx + (-v*dCy-C*dvy+B*dwy)*dtdy

              ! Cell-centered predicted states
              r = r + sr0
              u = u + su0
              v = v + sv0
              w = w + sw0
              p = p + sp0
              C = C + sC0
              A = 0.5*(AL+AR)
              B = 0.5*(BL+BR)

              ! Face averaged right state at left interface
              qp(l,i,j,k,ir,1) = r - drx
              qp(l,i,j,k,iu,1) = u - dux
              qp(l,i,j,k,iv,1) = v - dvx
              qp(l,i,j,k,iw,1) = w - dwx
              qp(l,i,j,k,ip,1) = p - dpx
              qp(l,i,j,k,iA,1) = AL     
              qp(l,i,j,k,iB,1) = B - dBx
              qp(l,i,j,k,iC,1) = C - dCx
              qp(l,i,j,k,ir,1) = MAX(smallr, qp(l,i,j,k,ir,1))
              qp(l,i,j,k,ip,1) = MAX(smallp, qp(l,i,j,k,ip,1))

              ! Face averaged left state at right interface
              qm(l,i,j,k,ir,1) = r + drx
              qm(l,i,j,k,iu,1) = u + dux
              qm(l,i,j,k,iv,1) = v + dvx
              qm(l,i,j,k,iw,1) = w + dwx
              qm(l,i,j,k,ip,1) = p + dpx
              qm(l,i,j,k,iA,1) = AR     
              qm(l,i,j,k,iB,1) = B + dBx
              qm(l,i,j,k,iC,1) = C + dCx
              qm(l,i,j,k,ir,1) = MAX(smallr, qm(l,i,j,k,ir,1))
              qm(l,i,j,k,ip,1) = MAX(smallp, qm(l,i,j,k,ip,1))

              ! Face averaged top state at bottom interface
              qp(l,i,j,k,ir,2) = r - dry
              qp(l,i,j,k,iu,2) = u - duy
              qp(l,i,j,k,iv,2) = v - dvy
              qp(l,i,j,k,iw,2) = w - dwy
              qp(l,i,j,k,ip,2) = p - dpy
              qp(l,i,j,k,iA,2) = A - dAy
              qp(l,i,j,k,iB,2) = BL     
              qp(l,i,j,k,iC,2) = C - dCy
              qp(l,i,j,k,ir,2) = MAX(smallr, qp(l,i,j,k,ir,2))
              qp(l,i,j,k,ip,2) = MAX(smallp, qp(l,i,j,k,ip,2))

              ! Face averaged bottom state at top interface
              qm(l,i,j,k,ir,2) = r + dry
              qm(l,i,j,k,iu,2) = u + duy
              qm(l,i,j,k,iv,2) = v + dvy
              qm(l,i,j,k,iw,2) = w + dwy
              qm(l,i,j,k,ip,2) = p + dpy
              qm(l,i,j,k,iA,2) = A + dAy
              qm(l,i,j,k,iB,2) = BR     
              qm(l,i,j,k,iC,2) = C + dCy
              qm(l,i,j,k,ir,2) = MAX(smallr, qm(l,i,j,k,ir,2))
              qm(l,i,j,k,ip,2) = MAX(smallp, qm(l,i,j,k,ip,2))

              ! Edge averaged right-top corner state (RT->LL)
              qRT(l,i,j,k,ir,3) = r + (+drx+dry)
              qRT(l,i,j,k,iu,3) = u + (+dux+duy)
              qRT(l,i,j,k,iv,3) = v + (+dvx+dvy)
              qRT(l,i,j,k,iw,3) = w + (+dwx+dwy)
              qRT(l,i,j,k,ip,3) = p + (+dpx+dpy)
              qRT(l,i,j,k,iC,3) = C + (+dCx+dCy)
              qRT(l,i,j,k,iA,3) = AR+ (   +dARy)
              qRT(l,i,j,k,iB,3) = BR+ (+dBRx   )
              qRT(l,i,j,k,ir,3) = MAX(smallr, qRT(l,i,j,k,ir,3))
              qRT(l,i,j,k,ip,3) = MAX(smallp, qRT(l,i,j,k,ip,3))

              ! Edge averaged right-bottom corner state (RB->LR)
              qRB(l,i,j,k,ir,3) = r + (+drx-dry)
              qRB(l,i,j,k,iu,3) = u + (+dux-duy)
              qRB(l,i,j,k,iv,3) = v + (+dvx-dvy)
              qRB(l,i,j,k,iw,3) = w + (+dwx-dwy)
              qRB(l,i,j,k,ip,3) = p + (+dpx-dpy)
              qRB(l,i,j,k,iC,3) = C + (+dCx-dCy)
              qRB(l,i,j,k,iA,3) = AR+ (   -dARy)
              qRB(l,i,j,k,iB,3) = BL+ (+dBLx   )
              qRB(l,i,j,k,ir,3) = MAX(smallr, qRB(l,i,j,k,ir,3))
              qRB(l,i,j,k,ip,3) = MAX(smallp, qRB(l,i,j,k,ip,3))

              ! Edge averaged left-top corner state (LT->RL)
              qLT(l,i,j,k,ir,3) = r + (-drx+dry)
              qLT(l,i,j,k,iu,3) = u + (-dux+duy)
              qLT(l,i,j,k,iv,3) = v + (-dvx+dvy)
              qLT(l,i,j,k,iw,3) = w + (-dwx+dwy)
              qLT(l,i,j,k,ip,3) = p + (-dpx+dpy)
              qLT(l,i,j,k,iC,3) = C + (-dCx+dCy)
              qLT(l,i,j,k,iA,3) = AL+ (   +dALy)
              qLT(l,i,j,k,iB,3) = BR+ (-dBRx   )
              qLT(l,i,j,k,ir,3) = MAX(smallr, qLT(l,i,j,k,ir,3))
              qLT(l,i,j,k,ip,3) = MAX(smallp, qLT(l,i,j,k,ip,3))

              ! Edge averaged left-bottom corner state (LB->RR)
              qLB(l,i,j,k,ir,3) = r + (-drx-dry)
              qLB(l,i,j,k,iu,3) = u + (-dux-duy)
              qLB(l,i,j,k,iv,3) = v + (-dvx-dvy)
              qLB(l,i,j,k,iw,3) = w + (-dwx-dwy)
              qLB(l,i,j,k,ip,3) = p + (-dpx-dpy)
              qLB(l,i,j,k,iC,3) = C + (-dCx-dCy)
              qLB(l,i,j,k,iA,3) = AL+ (   -dALy)
              qLB(l,i,j,k,iB,3) = BL+ (-dBLx   )
              qLB(l,i,j,k,ir,3) = MAX(smallr, qLB(l,i,j,k,ir,3))
              qLB(l,i,j,k,ip,3) = MAX(smallp, qLB(l,i,j,k,ip,3))

           END DO
        END DO
     END DO
  END DO

#if NVAR > 8
  ! Passive scalars
  DO n = 9, nvar
     DO k = klo, khi
        DO j = jlo, jhi
           DO i = ilo, ihi
              DO l = 1, ngrid
                 r   = q(l,i,j,k,n )              ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 v   = q(l,i,j,k,iv)
                 drx = half * dq(l,i,j,k,n,1)     ! TVD slopes
                 dry = half * dq(l,i,j,k,n,2)
                 sr0 = -u*drx*dtdx -v*dry*dtdy    ! Source terms
                 r   = r + sr0                    ! Predicted state
                 qp(l,i,j,k,n,1) = r - drx        ! Right state
                 qm(l,i,j,k,n,1) = r + drx        ! Left state
                 qp(l,i,j,k,n,2) = r - dry        ! Top state
                 qm(l,i,j,k,n,2) = r + dry        ! Bottom state
              END DO
           END DO
        END DO
     END DO
  END DO
#endif

END SUBROUTINE trace2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM==3
SUBROUTINE trace3d(q,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy,dz,dt,ngrid)
  USE amr_parameters
  USE hydro_parameters
  USE const
  IMPLICIT NONE

  INTEGER ::ngrid
  REAL(dp)::dx, dy, dz, dt

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  REAL(dp),DIMENSION(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:ndim)::dbf
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRB
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLB

  ! Declare local variables
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::Ex
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::Ey
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::Ez

  INTEGER ::i, j, k, l, n
  INTEGER ::ilo,ihi,jlo,jhi,klo,khi
  INTEGER ::ir, iu, iv, iw, ip, iA, iB, iC 
  REAL(dp)::dtdx, dtdy, dtdz, smallp
  REAL(dp)::r, u, v, w, p, A, B, C
  REAL(dp)::ELL, ELR, ERL, ERR
  REAL(dp)::FLL, FLR, FRL, FRR
  REAL(dp)::GLL, GLR, GRL, GRR
  REAL(dp)::drx, dux, dvx, dwx, dpx, dAx, dBx, dCx
  REAL(dp)::dry, duy, dvy, dwy, dpy, dAy, dbY, dCy
  REAL(dp)::drz, duz, dvz, dwz, dpz, dAz, dBz, dCz
  REAL(dp)::sr0, su0=0, sv0=0, sw0=0, sp0, sA0, sB0, sC0
  REAL(dp)::AL, AR, BL, BR, CL, CR
  REAL(dp)::dALy, dARy, dALz, dARz
  REAL(dp)::dBLx, dBRx, dBLz, dBRz
  REAL(dp)::dCLx, dCRx, dCLy, dCRy
  REAL(DP)::sAL0, sAR0, sBL0, sBR0, sCL0, sCR0
  
  dtdx = dt/dx
  dtdy = dt/dy
  dtdz = dt/dz
  smallp = smallr*smallc**2/gamma

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5; ia=6; ib=7; ic=8 

  DO k = klo, ku2
     DO j = jlo, ju2
        DO i = ilo, iu2
           DO l = 1, ngrid
              v = 0.25*(q(l,i,j-1,k-1,iv)+q(l,i,j-1,k,iv)+q(l,i,j,k-1,iv)+q(l,i,j,k,iv))
              w = 0.25*(q(l,i,j-1,k-1,iw)+q(l,i,j-1,k,iw)+q(l,i,j,k-1,iw)+q(l,i,j,k,iw))
              B = 0.5*(bf(l,i,j,k-1,2)+bf(l,i,j,k,2))
              C = 0.5*(bf(l,i,j-1,k,3)+bf(l,i,j,k,3))
              Ex(l,i,j,k)=v*C-w*B

              u = 0.25*(q(l,i-1,j,k-1,iu)+q(l,i-1,j,k,iu)+q(l,i,j,k-1,iu)+q(l,i,j,k,iu))
              w = 0.25*(q(l,i-1,j,k-1,iw)+q(l,i-1,j,k,iw)+q(l,i,j,k-1,iw)+q(l,i,j,k,iw))
              A = 0.5*(bf(l,i,j,k-1,1)+bf(l,i,j,k,1))
              C = 0.5*(bf(l,i-1,j,k,3)+bf(l,i,j,k,3))
              Ey(l,i,j,k)=w*A-u*C

              u = 0.25*(q(l,i-1,j-1,k,iu)+q(l,i-1,j,k,iu)+q(l,i,j-1,k,iu)+q(l,i,j,k,iu))
              v = 0.25*(q(l,i-1,j-1,k,iv)+q(l,i-1,j,k,iv)+q(l,i,j-1,k,iv)+q(l,i,j,k,iv))
              A = 0.5*(bf(l,i,j-1,k,1)+bf(l,i,j,k,1))
              B = 0.5*(bf(l,i-1,j,k,2)+bf(l,i,j,k,2))
              Ez(l,i,j,k)=u*B-v*A
           END DO
        END DO
     END DO
  END DO

  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi
           DO l = 1, ngrid

              ! Cell centered values
              r =    q(l,i,j,k,ir)
              u =    q(l,i,j,k,iu)
              v =    q(l,i,j,k,iv)
              w =    q(l,i,j,k,iw)            
              p =    q(l,i,j,k,ip)
              A =    q(l,i,j,k,iA)
              B =    q(l,i,j,k,iB)
              C =    q(l,i,j,k,iC)            

              ! Face centered variables
              AL =  bf(l,i  ,j  ,k  ,1)
              AR =  bf(l,i+1,j  ,k  ,1)
              BL =  bf(l,i  ,j  ,k  ,2)
              BR =  bf(l,i  ,j+1,k  ,2)
              CL =  bf(l,i  ,j  ,k  ,3)
              CR =  bf(l,i  ,j  ,k+1,3)

              ! Cell centered TVD slopes in X, Y and Z directions
              drx = half * dq(l,i,j,k,ir,1)
              dux = half * dq(l,i,j,k,iu,1)
              dvx = half * dq(l,i,j,k,iv,1)
              dwx = half * dq(l,i,j,k,iw,1)
              dpx = half * dq(l,i,j,k,ip,1)
              dBx = half * dq(l,i,j,k,iB,1)
              dCx = half * dq(l,i,j,k,iC,1)

              dry = half * dq(l,i,j,k,ir,2)
              duy = half * dq(l,i,j,k,iu,2)
              dvy = half * dq(l,i,j,k,iv,2)
              dwy = half * dq(l,i,j,k,iw,2)
              dpy = half * dq(l,i,j,k,ip,2)
              dAy = half * dq(l,i,j,k,iA,2)
              dCy = half * dq(l,i,j,k,iC,2)

              drz = half * dq(l,i,j,k,ir,3)
              duz = half * dq(l,i,j,k,iu,3)
              dvz = half * dq(l,i,j,k,iv,3)
              dwz = half * dq(l,i,j,k,iw,3)
              dpz = half * dq(l,i,j,k,ip,3)
              dAz = half * dq(l,i,j,k,iA,3)
              dBz = half * dq(l,i,j,k,iB,3)

              ! Face centered TVD slopes in transverse directions
              dALy = half * dbf(l,i  ,j  ,k  ,1,1)
              dARy = half * dbf(l,i+1,j  ,k  ,1,1)
              dALz = half * dbf(l,i  ,j  ,k  ,1,2)
              dARz = half * dbf(l,i+1,j  ,k  ,1,2)

              dBLx = half * dbf(l,i  ,j  ,k  ,2,1)
              dBRx = half * dbf(l,i  ,j+1,k  ,2,1)
              dBLz = half * dbf(l,i  ,j  ,k  ,2,2)
              dBRz = half * dbf(l,i  ,j+1,k  ,2,2)

              dCLx = half * dbf(l,i  ,j  ,k  ,3,1)
              dCRx = half * dbf(l,i  ,j  ,k+1,3,1)
              dCLy = half * dbf(l,i  ,j  ,k  ,3,2)
              dCRy = half * dbf(l,i  ,j  ,k+1,3,2)

              ! Edge centered electric field in X, Y and Z directions
              ELL = Ex(l,i,j  ,k  )
              ELR = Ex(l,i,j  ,k+1)
              ERL = Ex(l,i,j+1,k  )
              ERR = Ex(l,i,j+1,k+1)

              FLL = Ey(l,i  ,j,k  )
              FLR = Ey(l,i  ,j,k+1)
              FRL = Ey(l,i+1,j,k  )
              FRR = Ey(l,i+1,j,k+1)

              GLL = Ez(l,i  ,j  ,k)
              GLR = Ez(l,i  ,j+1,k)
              GRL = Ez(l,i+1,j  ,k)
              GRR = Ez(l,i+1,j+1,k)

              ! Face-centered predicted states
              sAL0 = +(GLR-GLL)*dtdy*half -(FLR-FLL)*dtdz*half
              sAR0 = +(GRR-GRL)*dtdy*half -(FRR-FRL)*dtdz*half
              sBL0 = -(GRL-GLL)*dtdx*half +(ELR-ELL)*dtdz*half
              sBR0 = -(GRR-GLR)*dtdx*half +(ERR-ERL)*dtdz*half
              sCL0 = +(FRL-FLL)*dtdx*half -(ERL-ELL)*dtdy*half
              sCR0 = +(FRR-FLR)*dtdx*half -(ERR-ELR)*dtdy*half
              
              AL = AL + sAL0
              AR = AR + sAR0
              BL = BL + sBL0
              BR = BR + sBR0
              CL = CL + sCL0
              CR = CR + sCR0
              
              ! Source terms (including transverse derivatives)
              sr0 = (-u*drx-dux*r)*dtdx + (-v*dry-dvy*r)*dtdy + (-w*drz-dwz*r)*dtdz 
              if(scheme.ne.'induction')then
              su0 = (-u*dux-(dpx+B*dBx+C*dCx)/r)*dtdx + (-v*duy+B*dAy/r)*dtdy + (-w*duz+C*dAz/r)*dtdz 
              sv0 = (-u*dvx+A*dBx/r)*dtdx + (-v*dvy-(dpy+A*dAy+C*dCy)/r)*dtdy + (-w*dvz+C*dBz/r)*dtdz
              sw0 = (-u*dwx+A*dCx/r)*dtdx + (-v*dwy+B*dCy/r)*dtdy + (-w*dwz-(dpz+A*dAz+B*dBz)/r)*dtdz 
              endif
              sp0 = (-u*dpx-dux*gamma*p)*dtdx + (-v*dpy-dvy*gamma*p)*dtdy + (-w*dpz-dwz*gamma*p)*dtdz

              ! Cell-centered predicted states
              r = r + sr0
              u = u + su0
              v = v + sv0
              w = w + sw0
              p = p + sp0
              A = 0.5*(AL+AR)
              B = 0.5*(BL+BR)
              C = 0.5*(CL+CR)

              ! Face averaged right state at left interface
              qp(l,i,j,k,ir,1) = r - drx
              qp(l,i,j,k,iu,1) = u - dux
              qp(l,i,j,k,iv,1) = v - dvx
              qp(l,i,j,k,iw,1) = w - dwx
              qp(l,i,j,k,ip,1) = p - dpx
              qp(l,i,j,k,iA,1) = AL
              qp(l,i,j,k,iB,1) = B - dBx
              qp(l,i,j,k,iC,1) = C - dCx
              qp(l,i,j,k,ir,1) = MAX(smallr, qp(l,i,j,k,ir,1))
              qp(l,i,j,k,ip,1) = MAX(smallp, qp(l,i,j,k,ip,1))

              ! Face averaged left state at right interface
              qm(l,i,j,k,ir,1) = r + drx
              qm(l,i,j,k,iu,1) = u + dux
              qm(l,i,j,k,iv,1) = v + dvx
              qm(l,i,j,k,iw,1) = w + dwx
              qm(l,i,j,k,ip,1) = p + dpx
              qm(l,i,j,k,iA,1) = AR
              qm(l,i,j,k,iB,1) = B + dBx
              qm(l,i,j,k,iC,1) = C + dCx
              qm(l,i,j,k,ir,1) = MAX(smallr, qm(l,i,j,k,ir,1))
              qm(l,i,j,k,ip,1) = MAX(smallp, qm(l,i,j,k,ip,1))

              ! Face averaged top state at bottom interface
              qp(l,i,j,k,ir,2) = r - dry
              qp(l,i,j,k,iu,2) = u - duy
              qp(l,i,j,k,iv,2) = v - dvy
              qp(l,i,j,k,iw,2) = w - dwy
              qp(l,i,j,k,ip,2) = p - dpy
              qp(l,i,j,k,iA,2) = A - dAy
              qp(l,i,j,k,iB,2) = BL
              qp(l,i,j,k,iC,2) = C - dCy
              qp(l,i,j,k,ir,2) = MAX(smallr, qp(l,i,j,k,ir,2))
              qp(l,i,j,k,ip,2) = MAX(smallp, qp(l,i,j,k,ip,2))

              ! Face averaged bottom state at top interface
              qm(l,i,j,k,ir,2) = r + dry
              qm(l,i,j,k,iu,2) = u + duy
              qm(l,i,j,k,iv,2) = v + dvy
              qm(l,i,j,k,iw,2) = w + dwy
              qm(l,i,j,k,ip,2) = p + dpy
              qm(l,i,j,k,iA,2) = A + dAy
              qm(l,i,j,k,iB,2) = BR
              qm(l,i,j,k,iC,2) = C + dCy
              qm(l,i,j,k,ir,2) = MAX(smallr, qm(l,i,j,k,ir,2))
              qm(l,i,j,k,ip,2) = MAX(smallp, qm(l,i,j,k,ip,2))

              ! Face averaged front state at back interface
              qp(l,i,j,k,ir,3) = r - drz
              qp(l,i,j,k,iu,3) = u - duz
              qp(l,i,j,k,iv,3) = v - dvz
              qp(l,i,j,k,iw,3) = w - dwz
              qp(l,i,j,k,ip,3) = p - dpz
              qp(l,i,j,k,iA,3) = A - dAz
              qp(l,i,j,k,iB,3) = B - dBz
              qp(l,i,j,k,iC,3) = CL
              qp(l,i,j,k,ir,3) = MAX(smallr, qp(l,i,j,k,ir,3))
              qp(l,i,j,k,ip,3) = MAX(smallp, qp(l,i,j,k,ip,3))

              ! Face averaged back state at front interface
              qm(l,i,j,k,ir,3) = r + drz
              qm(l,i,j,k,iu,3) = u + duz
              qm(l,i,j,k,iv,3) = v + dvz
              qm(l,i,j,k,iw,3) = w + dwz
              qm(l,i,j,k,ip,3) = p + dpz
              qm(l,i,j,k,iA,3) = A + dAz
              qm(l,i,j,k,iB,3) = B + dBz
              qm(l,i,j,k,iC,3) = CR
              qm(l,i,j,k,ir,3) = MAX(smallr, qm(l,i,j,k,ir,3))
              qm(l,i,j,k,ip,3) = MAX(smallp, qm(l,i,j,k,ip,3))

              ! X-edge averaged right-top corner state (RT->LL)
              qRT(l,i,j,k,ir,1) = r + (+dry+drz)
              qRT(l,i,j,k,iu,1) = u + (+duy+duz)
              qRT(l,i,j,k,iv,1) = v + (+dvy+dvz)
              qRT(l,i,j,k,iw,1) = w + (+dwy+dwz)
              qRT(l,i,j,k,ip,1) = p + (+dpy+dpz)
              qRT(l,i,j,k,iA,1) = A + (+dAy+dAz)
              qRT(l,i,j,k,iB,1) = BR+ (   +dBRz)
              qRT(l,i,j,k,iC,1) = CR+ (+dCRy   )
              qRT(l,i,j,k,ir,1) = MAX(smallr, qRT(l,i,j,k,ir,1))
              qRT(l,i,j,k,ip,1) = MAX(smallp, qRT(l,i,j,k,ip,1))

              ! X-edge averaged right-bottom corner state (RB->LR)
              qRB(l,i,j,k,ir,1) = r + (+dry-drz)
              qRB(l,i,j,k,iu,1) = u + (+duy-duz)
              qRB(l,i,j,k,iv,1) = v + (+dvy-dvz)
              qRB(l,i,j,k,iw,1) = w + (+dwy-dwz)
              qRB(l,i,j,k,ip,1) = p + (+dpy-dpz)
              qRB(l,i,j,k,iA,1) = A + (+dAy-dAz)
              qRB(l,i,j,k,iB,1) = BR+ (   -dBRz)
              qRB(l,i,j,k,iC,1) = CL+ (+dCLy   )
              qRB(l,i,j,k,ir,1) = MAX(smallr, qRB(l,i,j,k,ir,1))
              qRB(l,i,j,k,ip,1) = MAX(smallp, qRB(l,i,j,k,ip,1))

              ! X-edge averaged left-top corner state (LT->RL)
              qLT(l,i,j,k,ir,1) = r + (-dry+drz)
              qLT(l,i,j,k,iu,1) = u + (-duy+duz)
              qLT(l,i,j,k,iv,1) = v + (-dvy+dvz)
              qLT(l,i,j,k,iw,1) = w + (-dwy+dwz)
              qLT(l,i,j,k,ip,1) = p + (-dpy+dpz)
              qLT(l,i,j,k,iA,1) = A + (-dAy+dAz)
              qLT(l,i,j,k,iB,1) = BL+ (   +dBLz)
              qLT(l,i,j,k,iC,1) = CR+ (-dCRy   )
              qLT(l,i,j,k,ir,1) = MAX(smallr, qLT(l,i,j,k,ir,1))
              qLT(l,i,j,k,ip,1) = MAX(smallp, qLT(l,i,j,k,ip,1))

              ! X-edge averaged left-bottom corner state (LB->RR)
              qLB(l,i,j,k,ir,1) = r + (-dry-drz)
              qLB(l,i,j,k,iu,1) = u + (-duy-duz)
              qLB(l,i,j,k,iv,1) = v + (-dvy-dvz)
              qLB(l,i,j,k,iw,1) = w + (-dwy-dwz)
              qLB(l,i,j,k,ip,1) = p + (-dpy-dpz)
              qLB(l,i,j,k,iA,1) = A + (-dAy-dAz)
              qLB(l,i,j,k,iB,1) = BL+ (   -dBLz)
              qLB(l,i,j,k,iC,1) = CL+ (-dCLy   )
              qLB(l,i,j,k,ir,1) = MAX(smallr, qLB(l,i,j,k,ir,1))
              qLB(l,i,j,k,ip,1) = MAX(smallp, qLB(l,i,j,k,ip,1))

              ! Y-edge averaged right-top corner state (RT->LL)
              qRT(l,i,j,k,ir,2) = r + (+drx+drz)
              qRT(l,i,j,k,iu,2) = u + (+dux+duz)
              qRT(l,i,j,k,iv,2) = v + (+dvx+dvz)
              qRT(l,i,j,k,iw,2) = w + (+dwx+dwz)
              qRT(l,i,j,k,ip,2) = p + (+dpx+dpz)
              qRT(l,i,j,k,iA,2) = AR+ (   +dARz)
              qRT(l,i,j,k,iB,2) = B + (+dBx+dBz)
              qRT(l,i,j,k,iC,2) = CR+ (+dCRx   )
              qRT(l,i,j,k,ir,2) = MAX(smallr, qRT(l,i,j,k,ir,2))
              qRT(l,i,j,k,ip,2) = MAX(smallp, qRT(l,i,j,k,ip,2))

              ! Y-edge averaged right-bottom corner state (RB->LR)
              qRB(l,i,j,k,ir,2) = r + (+drx-drz)
              qRB(l,i,j,k,iu,2) = u + (+dux-duz)
              qRB(l,i,j,k,iv,2) = v + (+dvx-dvz)
              qRB(l,i,j,k,iw,2) = w + (+dwx-dwz)
              qRB(l,i,j,k,ip,2) = p + (+dpx-dpz)
              qRB(l,i,j,k,iA,2) = AR+ (   -dARz)
              qRB(l,i,j,k,iB,2) = B + (+dBx-dBz)
              qRB(l,i,j,k,iC,2) = CL+ (+dCLx   )
              qRB(l,i,j,k,ir,2) = MAX(smallr, qRB(l,i,j,k,ir,2))
              qRB(l,i,j,k,ip,2) = MAX(smallp, qRB(l,i,j,k,ip,2))

              ! Y-edge averaged left-top corner state (LT->RL)
              qLT(l,i,j,k,ir,2) = r + (-drx+drz)
              qLT(l,i,j,k,iu,2) = u + (-dux+duz)
              qLT(l,i,j,k,iv,2) = v + (-dvx+dvz)
              qLT(l,i,j,k,iw,2) = w + (-dwx+dwz)
              qLT(l,i,j,k,ip,2) = p + (-dpx+dpz)
              qLT(l,i,j,k,iA,2) = AL+ (   +dALz)
              qLT(l,i,j,k,iB,2) = B + (-dBx+dBz)
              qLT(l,i,j,k,iC,2) = CR+ (-dCRx   )
              qLT(l,i,j,k,ir,2) = MAX(smallr, qLT(l,i,j,k,ir,2))
              qLT(l,i,j,k,ip,2) = MAX(smallp, qLT(l,i,j,k,ip,2))

              ! Y-edge averaged left-bottom corner state (LB->RR)
              qLB(l,i,j,k,ir,2) = r + (-drx-drz)
              qLB(l,i,j,k,iu,2) = u + (-dux-duz)
              qLB(l,i,j,k,iv,2) = v + (-dvx-dvz)
              qLB(l,i,j,k,iw,2) = w + (-dwx-dwz)
              qLB(l,i,j,k,ip,2) = p + (-dpx-dpz)
              qLB(l,i,j,k,iA,2) = AL+ (   -dALz)
              qLB(l,i,j,k,iB,2) = B + (-dBx-dBz)
              qLB(l,i,j,k,iC,2) = CL+ (-dCLx   )
              qLB(l,i,j,k,ir,2) = MAX(smallr, qLB(l,i,j,k,ir,2))
              qLB(l,i,j,k,ip,2) = MAX(smallp, qLB(l,i,j,k,ip,2))

              ! Z-edge averaged right-top corner state (RT->LL)
              qRT(l,i,j,k,ir,3) = r + (+drx+dry)
              qRT(l,i,j,k,iu,3) = u + (+dux+duy)
              qRT(l,i,j,k,iv,3) = v + (+dvx+dvy)
              qRT(l,i,j,k,iw,3) = w + (+dwx+dwy)
              qRT(l,i,j,k,ip,3) = p + (+dpx+dpy)
              qRT(l,i,j,k,iA,3) = AR+ (   +dARy)
              qRT(l,i,j,k,iB,3) = BR+ (+dBRx   )
              qRT(l,i,j,k,iC,3) = C + (+dCx+dCy)
              qRT(l,i,j,k,ir,3) = MAX(smallr, qRT(l,i,j,k,ir,3))
              qRT(l,i,j,k,ip,3) = MAX(smallp, qRT(l,i,j,k,ip,3))

              ! Z-edge averaged right-bottom corner state (RB->LR)
              qRB(l,i,j,k,ir,3) = r + (+drx-dry)
              qRB(l,i,j,k,iu,3) = u + (+dux-duy)
              qRB(l,i,j,k,iv,3) = v + (+dvx-dvy)
              qRB(l,i,j,k,iw,3) = w + (+dwx-dwy)
              qRB(l,i,j,k,ip,3) = p + (+dpx-dpy)
              qRB(l,i,j,k,iA,3) = AR+ (   -dARy)
              qRB(l,i,j,k,iB,3) = BL+ (+dBLx   )
              qRB(l,i,j,k,iC,3) = C + (+dCx-dCy)
              qRB(l,i,j,k,ir,3) = MAX(smallr, qRB(l,i,j,k,ir,3))
              qRB(l,i,j,k,ip,3) = MAX(smallp, qRB(l,i,j,k,ip,3))

              ! Z-edge averaged left-top corner state (LT->RL)
              qLT(l,i,j,k,ir,3) = r + (-drx+dry)
              qLT(l,i,j,k,iu,3) = u + (-dux+duy)
              qLT(l,i,j,k,iv,3) = v + (-dvx+dvy)
              qLT(l,i,j,k,iw,3) = w + (-dwx+dwy)
              qLT(l,i,j,k,ip,3) = p + (-dpx+dpy)
              qLT(l,i,j,k,iA,3) = AL+ (   +dALy)
              qLT(l,i,j,k,iB,3) = BR+ (-dBRx   )
              qLT(l,i,j,k,iC,3) = C + (-dCx+dCy)
              qLT(l,i,j,k,ir,3) = MAX(smallr, qLT(l,i,j,k,ir,3))
              qLT(l,i,j,k,ip,3) = MAX(smallp, qLT(l,i,j,k,ip,3))

              ! Z-edge averaged left-bottom corner state (LB->RR)
              qLB(l,i,j,k,ir,3) = r + (-drx-dry)
              qLB(l,i,j,k,iu,3) = u + (-dux-duy)
              qLB(l,i,j,k,iv,3) = v + (-dvx-dvy)
              qLB(l,i,j,k,iw,3) = w + (-dwx-dwy)
              qLB(l,i,j,k,ip,3) = p + (-dpx-dpy)
              qLB(l,i,j,k,iA,3) = AL+ (   -dALy)
              qLB(l,i,j,k,iB,3) = BL+ (-dBLx   )
              qLB(l,i,j,k,iC,3) = C + (-dCx-dCy)
              qLB(l,i,j,k,ir,3) = MAX(smallr, qLB(l,i,j,k,ir,3))
              qLB(l,i,j,k,ip,3) = MAX(smallp, qLB(l,i,j,k,ip,3))

           END DO
        END DO
     END DO
  END DO

#if NVAR > 8
  ! Passive scalars
  DO n = 9, nvar
     DO k = klo, khi
        DO j = jlo, jhi
           DO i = ilo, ihi
              DO l = 1, ngrid
                 r   = q(l,i,j,k,n )            ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 v   = q(l,i,j,k,iv)
                 w   = q(l,i,j,k,iw)
                 drx = half * dq(l,i,j,k,n,1)   ! TVD slopes
                 dry = half * dq(l,i,j,k,n,2)
                 drz = half * dq(l,i,j,k,n,3)
                 sr0 = -u*drx*dtdx -v*dry*dtdy -w*drz*dtdz   ! Source terms
                 r   = r + sr0                  ! Predicted state
                 qp(l,i,j,k,n,1) = r - drx      ! Right state
                 qm(l,i,j,k,n,1) = r + drx      ! Left state
                 qp(l,i,j,k,n,2) = r - dry      ! Top state
                 qm(l,i,j,k,n,2) = r + dry      ! Bottom state
                 qp(l,i,j,k,n,3) = r - drz      ! Front state
                 qm(l,i,j,k,n,3) = r + drz      ! Back state
              END DO
           END DO
        END DO
     END DO
  END DO
#endif

END SUBROUTINE trace3d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpflxm(qm,im1,im2,jm1,jm2,km1,km2, &
     &             qp,ip1,ip2,jp1,jp2,kp1,kp2, &
     &                ilo,ihi,jlo,jhi,klo,khi, &
     &                ln ,lt1,lt2,bn ,bt1,bt2, flx,tmp,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  integer ::ln,lt1,lt2,bn,bt1,bt2
  integer ::im1,im2,jm1,jm2,km1,km2
  integer ::ip1,ip2,jp1,jp2,kp1,kp2
  integer ::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(1:nvector,im1:im2,jm1:jm2,km1:km2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar,1:ndim)::qp 
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar)::flx
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:2)::tmp
  
  ! local variables
  integer ::i, j, k, n, l, idim, xdim
  real(dp),dimension(1:nvar)::qleft,qright,fgdnv
  REAL(dp)::zero_flux, bn_mean, entho

  xdim=ln-1
  entho=one/(gamma-one)

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
           
              ! Enforce continuity for normal magnetic field
              bn_mean = half*(qm(l,i,j,k,bn,xdim)+qp(l,i,j,k,bn,xdim))

              ! Left state
              qleft (1) = qm(l,i,j,k,1  ,xdim) ! Mass density
              qleft (2) = qm(l,i,j,k,5  ,xdim) ! Pressure
              qleft (3) = qm(l,i,j,k,ln ,xdim) ! Normal velocity
              qleft (4) = bn_mean              ! Normal magnetic field
              qleft (5) = qm(l,i,j,k,lt1,xdim) ! Tangential velocity 1
              qleft (6) = qm(l,i,j,k,bt1,xdim) ! Tangential magnetic field 1
              qleft (7) = qm(l,i,j,k,lt2,xdim) ! Tangential velocity 2
              qleft (8) = qm(l,i,j,k,bt2,xdim) ! Tangential magnetic field 2

              ! Right state
              qright(1) = qp(l,i,j,k,1  ,xdim) ! Mass density
              qright(2) = qp(l,i,j,k,5  ,xdim) ! Pressure
              qright(3) = qp(l,i,j,k,ln ,xdim) ! Normal velocity
              qright(4) = bn_mean              ! Normal magnetic field
              qright(5) = qp(l,i,j,k,lt1,xdim) ! Tangential velocity 1
              qright(6) = qp(l,i,j,k,bt1,xdim) ! Tangential magnetic field 1
              qright(7) = qp(l,i,j,k,lt2,xdim) ! Tangential velocity 2
              qright(8) = qp(l,i,j,k,bt2,xdim) ! Tangential magnetic field 2

              ! Other advected quantities
#if NVAR > 8
              do n = 9, nvar
                 qleft (n) = qm(l,i,j,k,n,xdim)
                 qright(n) = qp(l,i,j,k,n,xdim)    
              end do
#endif
              ! Solve 1D Riemann problem
              zero_flux = one
              IF(scheme.NE.'induction')THEN
              SELECT CASE (riemann)
              CASE ('roe')
                 CALL athena_roe    (qleft,qright,fgdnv,zero_flux)
              CASE ('llf')
                 CALL lax_friedrich (qleft,qright,fgdnv,zero_flux)
              CASE ('hll')
                 CALL hll           (qleft,qright,fgdnv)
              CASE ('hlld')
                 CALL hlld          (qleft,qright,fgdnv)
              CASE ('upwind')
                 CALL lax_friedrich (qleft,qright,fgdnv,zero_flux)
              CASE ('hydro')
                 CALL hydro_acoustic(qleft,qright,fgdnv)
              CASE DEFAULT
                 write(*,*)'unknown riemann solver'
                 call clean_stop
              END SELECT
              ELSE
                 CALL upwind(qleft,qright,fgdnv,zero_flux)
              ENDIF
           
              ! Output fluxes
              flx(l,i,j,k,1  ) = fgdnv(1)  ! Mass density
              flx(l,i,j,k,5  ) = fgdnv(2)  ! Total energy
              flx(l,i,j,k,ln ) = fgdnv(3)  ! Normal momentum
              flx(l,i,j,k,bn ) = fgdnv(4)  ! Normal magnetic field
              flx(l,i,j,k,lt1) = fgdnv(5)  ! Transverse momentum 1
              flx(l,i,j,k,bt1) = fgdnv(6)  ! Transverse magnetic field 1
              flx(l,i,j,k,lt2) = fgdnv(7)  ! Transverse momentum 2
              flx(l,i,j,k,bt2) = fgdnv(8)  ! Transverse magnetic field 2

              ! Other advected quantities
#if NVAR > 8
              do n = 9, nvar
                 flx(l,i,j,k,n) = fgdnv(n)
              end do
#endif  
              ! Normal velocity estimate
              tmp(l,i,j,k,1) = half*(qleft(3)+qright(3))
              ! Internal energy flux
              if(fgdnv(1)>zero)then
                 tmp(l,i,j,k,2) = qleft (2)/qleft (1)*entho*fgdnv(1)
              else
                 tmp(l,i,j,k,2) = qright(2)/qright(1)*entho*fgdnv(1)
              end if

           end do
        end do
     end do
  end do
  
end subroutine cmpflxm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE cmp_mag_flx(qRT,irt1,irt2,jrt1,jrt2,krt1,krt2, &
       &               qRB,irb1,irb2,jrb1,jrb2,krb1,krb2, &
       &               qLT,ilt1,ilt2,jlt1,jlt2,klt1,klt2, &
       &               qLB,ilb1,ilb2,jlb1,jlb2,klb1,klb2, &
       &                   ilo ,ihi ,jlo ,jhi ,klo ,khi , &
       &                   lp1 ,lp2 ,lor ,bp1 ,bp2 ,bor ,emf,ngrid)
  ! 2D Riemann solver to compute EMF at cell edges
  USE amr_parameters
  USE hydro_parameters
  USE const
  IMPLICIT NONE

  INTEGER ::ngrid
  ! indices of the 2 planar velocity lp1 and lp2 and the orthogonal one,
  ! lor and idem for the magnetic field
  INTEGER ::lp1,lp2,lor,bp1,bp2,bor
  INTEGER ::irt1,irt2,jrt1,jrt2,krt1,krt2
  INTEGER ::irb1,irb2,jrb1,jrb2,krb1,krb2
  INTEGER ::ilt1,ilt2,jlt1,jlt2,klt1,klt2
  INTEGER ::ilb1,ilb2,jlb1,jlb2,klb1,klb2
  INTEGER ::ilo,ihi,jlo,jhi,klo,khi
  REAL(dp),DIMENSION(1:nvector,irt1:irt2,jrt1:jrt2,krt1:krt2,1:nvar,1:3)::qRT
  REAL(dp),DIMENSION(1:nvector,irb1:irb2,jrb1:jrb2,krb1:krb2,1:nvar,1:3)::qRB
  REAL(dp),DIMENSION(1:nvector,ilt1:ilt2,jlt1:jlt2,klt1:klt2,1:nvar,1:3)::qLT
  REAL(dp),DIMENSION(1:nvector,ilb1:ilb2,jlb1:jlb2,klb1:klb2,1:nvar,1:3)::qLB

  REAL(dp),DIMENSION(1:nvector,ilb1:ilb2,jlb1:jlb2,klb1:klb2):: emf

  ! local variables
  INTEGER ::i, j, k, n, l, idim, xdim, m
  REAL(dp),DIMENSION(1:nvector,1:nvar)::qLL,qRL,qLR,qRR
  REAL(dp),DIMENSION(1:nvar)::qleft,qright,fmean_x,fmean_y,qtmp
  REAL(dp) :: ELL,ERL,ELR,ERR,SL,SR,SB,ST,SAL,SAR,SAT,SAB
  REAL(dp) :: zero_flux,E
  REAL(dp) :: cLLx,cRLx,cLRx,cRRx,cLLy,cRLy,cLRy,cRRy
  REAL(dp) :: cfastLLx,cfastRLx,cfastLRx,cfastRRx,cfastLLy,cfastRLy,cfastLRy,cfastRRy
  REAL(dp) :: calfvenR,calfvenL,calfvenT,calfvenB
  REAL(dp) :: vLLx,vRLx,vLRx,vRRx,vLLy,vRLy,vLRy,vRRy
  REAL(dp) :: rLL,rLR,rRL,rRR,pLL,pLR,pRL,pRR,uLL,uLR,uRL,uRR,vLL,vLR,vRL,vRR
  REAL(dp) :: ALL,ALR,ARL,ARR,BLL,BLR,BRL,BRR,CLL,CLR,CRL,CRR
  REAL(dp) :: PtotLL,PtotLR,PtotRL,PtotRR,rcLLx,rcLRx,rcRLx,rcRRx,rcLLy,rcLRy,rcRLy,rcRRy
  REAL(dp) :: ustar,vstar,rstarLLx,rstarLRx,rstarRLx,rstarRRx,rstarLLy,rstarLRy,rstarRLy,rstarRRy
  REAL(dp) :: rstarLL,rstarLR,rstarRL,rstarRR,AstarLL,AstarLR,AstarRL,AstarRR,BstarLL,BstarLR,BstarRL,BstarRR
  REAL(dp) :: EstarLLx,EstarLRx,EstarRLx,EstarRRx,EstarLLy,EstarLRy,EstarRLy,EstarRRy,EstarLL,EstarLR,EstarRL,EstarRR
  REAL(dp) :: AstarT,AstarB,BstarR,BstarL

  xdim = lor - 1

  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi

           ! Density
           DO l = 1, ngrid
              qLL (l,1) = qRT(l,i,j,k,1,xdim)
              qRL (l,1) = qLT(l,i,j,k,1,xdim)
              qLR (l,1) = qRB(l,i,j,k,1,xdim)
              qRR (l,1) = qLB(l,i,j,k,1,xdim)
           END DO           

           ! Pressure 
           DO l = 1, ngrid
              qLL (l,2) = qRT(l,i,j,k,5,xdim)
              qRL (l,2) = qLT(l,i,j,k,5,xdim)
              qLR (l,2) = qRB(l,i,j,k,5,xdim)
              qRR (l,2) = qLB(l,i,j,k,5,xdim)
           END DO

           ! First parallel velocity 
           DO l = 1, ngrid
              qLL (l,3) = qRT(l,i,j,k,lp1,xdim)
              qRL (l,3) = qLT(l,i,j,k,lp1,xdim)
              qLR (l,3) = qRB(l,i,j,k,lp1,xdim)
              qRR (l,3) = qLB(l,i,j,k,lp1,xdim)
           END DO

           ! Second parallel velocity 
           DO l = 1, ngrid
              qLL (l,4) = qRT(l,i,j,k,lp2,xdim)
              qRL (l,4) = qLT(l,i,j,k,lp2,xdim)
              qLR (l,4) = qRB(l,i,j,k,lp2,xdim)
              qRR (l,4) = qLB(l,i,j,k,lp2,xdim)
           END DO

           ! First parallel magnetic field (enforce continuity)
           DO l = 1, ngrid
              qLL (l,6) = half*(qRT(l,i,j,k,bp1,xdim)+qLT(l,i,j,k,bp1,xdim))
              qRL (l,6) = half*(qRT(l,i,j,k,bp1,xdim)+qLT(l,i,j,k,bp1,xdim))
              qLR (l,6) = half*(qRB(l,i,j,k,bp1,xdim)+qLB(l,i,j,k,bp1,xdim))
              qRR (l,6) = half*(qRB(l,i,j,k,bp1,xdim)+qLB(l,i,j,k,bp1,xdim))
           END DO

           ! Second parallel magnetic field (enforce continuity)
           DO l = 1, ngrid
              qLL (l,7) = half*(qRT(l,i,j,k,bp2,xdim)+qRB(l,i,j,k,bp2,xdim))
              qRL (l,7) = half*(qLT(l,i,j,k,bp2,xdim)+qLB(l,i,j,k,bp2,xdim))
              qLR (l,7) = half*(qRT(l,i,j,k,bp2,xdim)+qRB(l,i,j,k,bp2,xdim))
              qRR (l,7) = half*(qLT(l,i,j,k,bp2,xdim)+qLB(l,i,j,k,bp2,xdim))
           END DO

           ! Orthogonal velocity 
           DO l = 1, ngrid
              qLL (l,5) = qRT(l,i,j,k,lor,xdim)
              qRL (l,5) = qLT(l,i,j,k,lor,xdim)
              qLR (l,5) = qRB(l,i,j,k,lor,xdim)
              qRR (l,5) = qLB(l,i,j,k,lor,xdim)
           END DO

           ! Orthogonal magnetic Field
           DO l = 1, ngrid
              qLL (l,8) = qRT(l,i,j,k,bor,xdim)
              qRL (l,8) = qLT(l,i,j,k,bor,xdim)
              qLR (l,8) = qRB(l,i,j,k,bor,xdim)
              qRR (l,8) = qLB(l,i,j,k,bor,xdim)
           END DO

           ! Compute final fluxes
            DO l = 1, ngrid 

               ! vx*by - vy*bx at the four edge centers
               ELL = qLL(l,3)*qLL(l,7) - qLL(l,4)*qLL(l,6)
               ERL = qRL(l,3)*qRL(l,7) - qRL(l,4)*qRL(l,6)
               ELR = qLR(l,3)*qLR(l,7) - qLR(l,4)*qLR(l,6)  
               ERR = qRR(l,3)*qRR(l,7) - qRR(l,4)*qRR(l,6) 

               if(riemann2d=='hlld')then

                  
                  rLL=qLL(l,1); pLL=qLL(l,2); uLL=qLL(l,3); vLL=qLL(l,4); ALL=qLL(l,6); BLL=qLL(l,7) ; CLL=qLL(l,8) 
                  rLR=qLR(l,1); pLR=qLR(l,2); uLR=qLR(l,3); vLR=qLR(l,4); ALR=qLR(l,6); BLR=qLR(l,7) ; CLR=qLR(l,8) 
                  rRL=qRL(l,1); pRL=qRL(l,2); uRL=qRL(l,3); vRL=qRL(l,4); ARL=qRL(l,6); BRL=qRL(l,7) ; CRL=qRL(l,8) 
                  rRR=qRR(l,1); pRR=qRR(l,2); uRR=qRR(l,3); vRR=qRR(l,4); ARR=qRR(l,6); BRR=qRR(l,7) ; CRR=qRR(l,8) 

                  ! Compute 4 fast magnetosonic velocity relative to x direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,3); qtmp(4)=qLL(l,6); qtmp(5)=qLL(l,4); qtmp(6)=qLL(l,7)
                  call find_speed_fast(qtmp,cfastLLx)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,3); qtmp(4)=qLR(l,6); qtmp(5)=qLR(l,4); qtmp(6)=qLR(l,7)
                  call find_speed_fast(qtmp,cfastLRx)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,3); qtmp(4)=qRL(l,6); qtmp(5)=qRL(l,4); qtmp(6)=qRL(l,7)
                  call find_speed_fast(qtmp,cfastRLx)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,3); qtmp(4)=qRR(l,6); qtmp(5)=qRR(l,4); qtmp(6)=qRR(l,7)
                  call find_speed_fast(qtmp,cfastRRx)

                  ! Compute 4 fast magnetosonic velocity relative to y direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,4); qtmp(4)=qLL(l,7); qtmp(5)=qLL(l,3); qtmp(6)=qLL(l,6)
                  call find_speed_fast(qtmp,cfastLLy)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,4); qtmp(4)=qLR(l,7); qtmp(5)=qLR(l,3); qtmp(6)=qLR(l,6)
                  call find_speed_fast(qtmp,cfastLRy)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,4); qtmp(4)=qRL(l,7); qtmp(5)=qRL(l,3); qtmp(6)=qRL(l,6)
                  call find_speed_fast(qtmp,cfastRLy)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,4); qtmp(4)=qRR(l,7); qtmp(5)=qRR(l,3); qtmp(6)=qRR(l,6)
                  call find_speed_fast(qtmp,cfastRRy)

                  SL=min(uLL,uLR,uRL,uRR)-max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
                  SR=max(uLL,uLR,uRL,uRR)+max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
                  SB=min(vLL,vLR,vRL,vRR)-max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
                  ST=max(vLL,vLR,vRL,vRR)+max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)

                  ELL=uLL*BLL-vLL*ALL
                  ELR=uLR*BLR-vLR*ALR
                  ERL=uRL*BRL-vRL*ARL
                  ERR=uRR*BRR-vRR*ARR

                  PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
                  PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
                  PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
                  PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
                  
                  rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL) 
                  rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
                  rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR) 
                  rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
                  
                  ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
                  vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)
                  
                  rstarLLx=rLL*(SL-uLL)/(SL-ustar); BstarLL=BLL*(SL-uLL)/(SL-ustar)
                  rstarLLy=rLL*(SB-vLL)/(SB-vstar); AstarLL=ALL*(SB-vLL)/(SB-vstar)
                  rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
                  EstarLLx=ustar*BstarLL-vLL  *ALL
                  EstarLLy=uLL  *BLL    -vstar*AstarLL
                  EstarLL =ustar*BstarLL-vstar*AstarLL

                  rstarLRx=rLR*(SL-uLR)/(SL-ustar); BstarLR=BLR*(SL-uLR)/(SL-ustar)
                  rstarLRy=rLR*(ST-vLR)/(ST-vstar); AstarLR=ALR*(ST-vLR)/(ST-vstar)
                  rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
                  EstarLRx=ustar*BstarLR-vLR  *ALR
                  EstarLRy=uLR  *BLR    -vstar*AstarLR
                  EstarLR =ustar*BstarLR-vstar*AstarLR

                  rstarRLx=rRL*(SR-uRL)/(SR-ustar); BstarRL=BRL*(SR-uRL)/(SR-ustar)
                  rstarRLy=rRL*(SB-vRL)/(SB-vstar); AstarRL=ARL*(SB-vRL)/(SB-vstar)
                  rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
                  EstarRLx=ustar*BstarRL-vRL  *ARL
                  EstarRLy=uRL  *BRL    -vstar*AstarRL
                  EstarRL =ustar*BstarRL-vstar*AstarRL

                  rstarRRx=rRR*(SR-uRR)/(SR-ustar); BstarRR=BRR*(SR-uRR)/(SR-ustar)
                  rstarRRy=rRR*(ST-vRR)/(ST-vstar); AstarRR=ARR*(ST-vRR)/(ST-vstar)
                  rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
                  EstarRRx=ustar*BstarRR-vRR  *ARR
                  EstarRRy=uRR  *BRR    -vstar*AstarRR
                  EstarRR =ustar*BstarRR-vstar*AstarRR

                  calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                       &       abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
                  calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                       &       abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)
                  calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                       &       abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
                  calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                       &       abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)
                  SAL=min(ustar-calfvenL,zero); SAR=max(ustar+calfvenR,zero)
                  SAB=min(vstar-calfvenB,zero); SAT=max(vstar+calfvenT,zero)
                  AstarT=(SAR*AstarRR-SAL*AstarLR)/(SAR-SAL); AstarB=(SAR*AstarRL-SAL*AstarLL)/(SAR-SAL)
                  BstarR=(SAT*BstarRR-SAB*BstarRL)/(SAT-SAB); BstarL=(SAT*BstarLR-SAB*BstarLL)/(SAT-SAB)

                  if(SB>0d0)then
                     if(SL>0d0)then
                     E=ELL
                     else if(SR<0d0)then
                     E=ERL
                     else
                     E=(SAR*EstarLLx-SAL*EstarRLx+SAR*SAL*(BRL-BLL))/(SAR-SAL)
                     endif
                  else if (ST<0d0)then
                     if(SL>0d0)then
                     E=ELR
                     else if(SR<0d0)then
                     E=ERR
                     else
                     E=(SAR*EstarLRx-SAL*EstarRRx+SAR*SAL*(BRR-BLR))/(SAR-SAL)
                     endif
                  else if(SL>0d0)then
                     E=(SAT*EstarLLy-SAB*EstarLRy-SAT*SAB*(ALR-ALL))/(SAT-SAB)
                  else if(SR<0d0)then
                     E=(SAT*EstarRLy-SAB*EstarRRy-SAT*SAB*(ARR-ARL))/(SAT-SAB)
                  else
                     E=(SAL*SAB*EstarRR-SAL*SAT*EstarRL-SAR*SAB*EstarLR+SAR*SAT*EstarLL)/(SAR-SAL)/(SAT-SAB) &
                          & -SAT*SAB/(SAT-SAB)*(AstarT-AstarB)+SAR*SAL/(SAR-SAL)*(BstarR-BstarL)
                  endif
                  
                  emf(l,i,j,k) = E

               else if(riemann2d=='hll')then

                  ! Compute 4 fast magnetosonic velocity relative to x direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,3); qtmp(4)=qLL(l,6); qtmp(5)=qLL(l,4); qtmp(6)=qLL(l,7)
                  vLLx=qtmp(3); call find_speed_fast(qtmp,cLLx)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,3); qtmp(4)=qLR(l,6); qtmp(5)=qLR(l,4); qtmp(6)=qLR(l,7)
                  vLRx=qtmp(3); call find_speed_fast(qtmp,cLRx)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,3); qtmp(4)=qRL(l,6); qtmp(5)=qRL(l,4); qtmp(6)=qRL(l,7)
                  vRLx=qtmp(3); call find_speed_fast(qtmp,cRLx)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,3); qtmp(4)=qRR(l,6); qtmp(5)=qRR(l,4); qtmp(6)=qRR(l,7)
                  vRRx=qtmp(3); call find_speed_fast(qtmp,cRRx)

                  ! Compute 4 fast magnetosonic velocity relative to y direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,4); qtmp(4)=qLL(l,7); qtmp(5)=qLL(l,3); qtmp(6)=qLL(l,6)
                  vLLy=qtmp(3); call find_speed_fast(qtmp,cLLy)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,4); qtmp(4)=qLR(l,7); qtmp(5)=qLR(l,3); qtmp(6)=qLR(l,6)
                  vLRy=qtmp(3); call find_speed_fast(qtmp,cLRy)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,4); qtmp(4)=qRL(l,7); qtmp(5)=qRL(l,3); qtmp(6)=qRL(l,6)
                  vRLy=qtmp(3); call find_speed_fast(qtmp,cRLy)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,4); qtmp(4)=qRR(l,7); qtmp(5)=qRR(l,3); qtmp(6)=qRR(l,6)
                  vRRy=qtmp(3); call find_speed_fast(qtmp,cRRy)

                  SL=min(min(vLLx,vLRx,VRLx,vRRx)-max(cLLx,cLRx,cRLx,cRRx),zero)
                  SR=max(max(vLLx,vLRx,VRLx,vRRx)+max(cLLx,cLRx,cRLx,cRRx),zero)
                  SB=min(min(vLLy,vLRy,VRLy,vRRy)-max(cLLy,cLRy,cRLy,cRRy),zero)
                  ST=max(max(vLLy,vLRy,VRLy,vRRy)+max(cLLy,cLRy,cRLy,cRRy),zero)

                  emf(l,i,j,k) = (SL*SB*ERR-SL*ST*ERL-SR*SB*ELR+SR*ST*ELL)/(SR-SL)/(ST-SB) &
                       -ST*SB/(ST-SB)*(qRR(l,6)-qLL(l,6)) &
                       +SR*SL/(SR-SL)*(qRR(l,7)-qLL(l,7))

               else if (riemann2d=='hlla')then

                  ! Compute 4 Alfven velocity relative to x direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,3); qtmp(4)=qLL(l,6); qtmp(5)=qLL(l,4); qtmp(6)=qLL(l,7)
                  vLLx=qtmp(3); call find_speed_alfven(qtmp,cLLx)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,3); qtmp(4)=qLR(l,6); qtmp(5)=qLR(l,4); qtmp(6)=qLR(l,7)
                  vLRx=qtmp(3); call find_speed_alfven(qtmp,cLRx)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,3); qtmp(4)=qRL(l,6); qtmp(5)=qRL(l,4); qtmp(6)=qRL(l,7)
                  vRLx=qtmp(3); call find_speed_alfven(qtmp,cRLx)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,3); qtmp(4)=qRR(l,6); qtmp(5)=qRR(l,4); qtmp(6)=qRR(l,7)
                  vRRx=qtmp(3); call find_speed_alfven(qtmp,cRRx)

                  ! Compute 4 Alfven relative to y direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,4); qtmp(4)=qLL(l,7); qtmp(5)=qLL(l,3); qtmp(6)=qLL(l,6)
                  vLLy=qtmp(3); call find_speed_alfven(qtmp,cLLy)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,4); qtmp(4)=qLR(l,7); qtmp(5)=qLR(l,3); qtmp(6)=qLR(l,6)
                  vLRy=qtmp(3); call find_speed_alfven(qtmp,cLRy)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,4); qtmp(4)=qRL(l,7); qtmp(5)=qRL(l,3); qtmp(6)=qRL(l,6)
                  vRLy=qtmp(3); call find_speed_alfven(qtmp,cRLy)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,4); qtmp(4)=qRR(l,7); qtmp(5)=qRR(l,3); qtmp(6)=qRR(l,6)
                  vRRy=qtmp(3); call find_speed_alfven(qtmp,cRRy)

                  SL=min(min(vLLx,vLRx,VRLx,vRRx)-max(cLLx,cLRx,cRLx,cRRx),zero)
                  SR=max(max(vLLx,vLRx,VRLx,vRRx)+max(cLLx,cLRx,cRLx,cRRx),zero)
                  SB=min(min(vLLy,vLRy,VRLy,vRRy)-max(cLLy,cLRy,cRLy,cRRy),zero)
                  ST=max(max(vLLy,vLRy,VRLy,vRRy)+max(cLLy,cLRy,cRLy,cRRy),zero)

                  emf(l,i,j,k) = (SL*SB*ERR-SL*ST*ERL-SR*SB*ELR+SR*ST*ELL)/(SR-SL)/(ST-SB) &
                       -ST*SB/(ST-SB)*(qRR(l,6)-qLL(l,6)) &
                       +SR*SL/(SR-SL)*(qRR(l,7)-qLL(l,7))

               else

                  ! find the average value of E
                  E = forth*(ELL+ERL+ELR+ERR)
                  
                  ! call the first solver in the x direction
                  ! density
                  qleft (1) = half*(qLL(l,1)+qLR(l,1))
                  qright(1) = half*(qRR(l,1)+qRL(l,1))

                  ! pressure
                  qleft (2) = half*(qLL(l,2)+qLR(l,2))
                  qright(2) = half*(qRR(l,2)+qRL(l,2))
                  
                  ! vt1 becomes normal velocity
                  qleft (3) = half*(qLL(l,3)+qLR(l,3))
                  qright(3) = half*(qRR(l,3)+qRL(l,3))
                  
                  ! bt1 becomes normal magnetic field
                  qleft (4) = half*(qLL(l,6)+qLR(l,6))
                  qright(4) = half*(qRR(l,6)+qRL(l,6))
                  
                  ! vt2 becomes transverse velocity field
                  qleft (5) = half*(qLL(l,4)+qLR(l,4))
                  qright(5) = half*(qRR(l,4)+qRL(l,4))
                  
                  ! bt2 becomes transverse magnetic field 
                  qleft (6) = half*(qLL(l,7)+qLR(l,7))
                  qright(6) = half*(qRR(l,7)+qRL(l,7))
                  
                  ! velocity component perp. to the plane is now transverse
                  qleft (7) = half*(qLL(l,5)+qLR(l,5))
                  qright(7) = half*(qRR(l,5)+qRL(l,5))
                  
                  ! magnetic field component perp. to the plane is now transverse
                  qleft (8) = half*(qLL(l,8)+qLR(l,8))
                  qright(8) = half*(qRR(l,8)+qRL(l,8))
                  
                  zero_flux = 0.0
                  SELECT CASE (riemann2d)
                  CASE ('roe')
                     CALL athena_roe   (qleft,qright,fmean_x,zero_flux)
                  CASE ('llf')
                     CALL lax_friedrich(qleft,qright,fmean_x,zero_flux)
                  CASE ('upwind')
                     CALL upwind       (qleft,qright,fmean_x,zero_flux)
                  CASE DEFAULT
                     write(*,*)'unknown 2D riemann solver'
                     call clean_stop
                  END SELECT

                  ! call the second solver in the y direction
                  ! density
                  qleft (1) = half*(qLL(l,1)+qRL(l,1))
                  qright(1) = half*(qRR(l,1)+qLR(l,1))
                  
                  ! pressure
                  qleft (2) = half*(qLL(l,2)+qRL(l,2))
                  qright(2) = half*(qRR(l,2)+qLR(l,2))
                  
                  ! vt2 becomes normal velocity
                  qleft (3) = half*(qLL(l,4)+qRL(l,4))
                  qright(3) = half*(qRR(l,4)+qLR(l,4))
                  
                  ! bt2 becomes normal magnetic field
                  qleft (4) = half*(qLL(l,7)+qRL(l,7))
                  qright(4) = half*(qRR(l,7)+qLR(l,7))
                  
                  ! vt1 becomes transverse velocity field 
                  qleft (5) = half*(qLL(l,3)+qRL(l,3))
                  qright(5) = half*(qRR(l,3)+qLR(l,3))
                  
                  ! bt1 becomes transverse magnetic field 
                  qleft (6) = half*(qLL(l,6)+qRL(l,6))
                  qright(6) = half*(qRR(l,6)+qLR(l,6))
                  
                  ! velocity component perp. to the plane is now transverse
                  qleft (7) = half*(qLL(l,5)+qRL(l,5))
                  qright(7) = half*(qRR(l,5)+qLR(l,5))
                  
                  ! magnetic field component perp. to the plane is now transverse
                  qleft (8) = half*(qLL(l,8)+qRL(l,8))
                  qright(8) = half*(qRR(l,8)+qLR(l,8))
                  
                  zero_flux = 0.
                  SELECT CASE (riemann2d)
                  CASE ('roe')
                     CALL athena_roe   (qleft,qright,fmean_y,zero_flux)
                  CASE ('llf')
                     CALL lax_friedrich(qleft,qright,fmean_y,zero_flux)
                  CASE ('upwind')
                     CALL upwind       (qleft,qright,fmean_y,zero_flux)
                  CASE DEFAULT
                     write(*,*)'unknown 2D riemann solver'
                     call clean_stop
                  END SELECT
                  
                  ! compute the final value of E including the 2D diffusive
                  ! terms that ensure stability 
                  emf(l,i,j,k) = E + (fmean_x(6) - fmean_y(6))
                  
               endif

            END DO
        END DO
     END DO
  END DO


END SUBROUTINE cmp_mag_flx
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim(uin,q,bf,gravin,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf  

  integer ::i, j, k, l, n, idim
  real(dp)::eint, smalle, smallp, etot
  real(dp),dimension(1:nvector),save::eken,emag

  smalle = smallc**2/gamma/(gamma-one)
  smallp = smallr*smallc**2/gamma

  ! Store face centered magnetic field
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2+1
           DO l = 1, ngrid
              if(i<=iu2)then
                 bf(l,i,j,k,1) = uin(l,i,j,k,6)
              else
                 bf(l,i,j,k,1) = uin(l,i-1,j,k,nvar+1)
              endif
           END DO
        end do
     end do
  end do
  do k = ku1, ku2
     do j = ju1, ju2+1
        do i = iu1, iu2
           DO l = 1, ngrid
              if(j<=ju2)then
                 bf(l,i,j,k,2) = uin(l,i,j,k,7)
              else
                 bf(l,i,j,k,2) = uin(l,i,j-1,k,nvar+2)
              endif
           END DO
        end do
     end do
  end do
  do k = ku1, ku2+1
     do j = ju1, ju2
        do i = iu1, iu2
           DO l = 1, ngrid
              if(k<=ku2)then
                 bf(l,i,j,k,3) = uin(l,i,j,k,8)
              else
                 bf(l,i,j,k,3) = uin(l,i,j,k-1,nvar+3)
              endif
           END DO
        end do
     end do
  end do

  ! Convert to primitive variable
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2

           ! Compute density
           do l = 1, ngrid
              q(l,i,j,k,1) = max(uin(l,i,j,k,1),smallr)
           end do

           ! Debug
           if(debug)then
              do l = 1, ngrid
                 if(uin(l,i,j,k,1).le.smallr)then
                    write(*,*)'negative density'
                    write(*,*)uin(l,i,j,k,1)
                    stop
                 end if
              end do
           end if

           ! Compute velocities
           do l = 1, ngrid
              q(l,i,j,k,2) = uin(l,i,j,k,2)/uin(l,i,j,k,1)
              q(l,i,j,k,3) = uin(l,i,j,k,3)/uin(l,i,j,k,1)
              q(l,i,j,k,4) = uin(l,i,j,k,4)/uin(l,i,j,k,1)
           end do

           ! Compute cell centered magnetic field
           DO l = 1, ngrid
              q(l,i,j,k,6) = (uin(l,i,j,k,6)+uin(l,i,j,k,nvar+1))*half
              q(l,i,j,k,7) = (uin(l,i,j,k,7)+uin(l,i,j,k,nvar+2))*half
              q(l,i,j,k,8) = (uin(l,i,j,k,8)+uin(l,i,j,k,nvar+3))*half
           END DO

           ! Compute specific kinetic energy and magnetic energy
           do l = 1, ngrid
              eken(l) = half*(q(l,i,j,k,2)**2+q(l,i,j,k,3)**2+q(l,i,j,k,4)**2)
              emag(l) = half*(q(l,i,j,k,6)**2+q(l,i,j,k,7)**2+q(l,i,j,k,8)**2)
           end do

           ! Compute thermal pressure through EOS
           do l = 1, ngrid
              etot = uin(l,i,j,k,5) - emag(l)
              eint = etot/uin(l,i,j,k,1)-eken(l)
              q(l,i,j,k,5)=MAX((gamma-one)*q(l,i,j,k,1)*eint,smallp)
           end do

           ! Gravity predictor step
           do idim = 1, ndim
              do l = 1, ngrid
                 q(l,i,j,k,idim+1) = q(l,i,j,k,idim+1) + gravin(l,i,j,k,idim)*dt*half
              end do
           end do

        end do
     end do
  end do

  ! Passive scalar
#if NVAR > 8
  do n = 9, nvar
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              do l = 1, ngrid
                 q(l,i,j,k,n) = uin(l,i,j,k,n)/uin(l,i,j,k,1)
              end do
           end do
        end do
     end do
  end do
#endif
 
end subroutine ctoprim
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine uslope(bf,q,dq,dbf,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer::ngrid
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:ndim)::dbf

  ! local arrays
  integer::i, j, k, l, n
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff,xslope_type
  integer::ilo,ihi,jlo,jhi,klo,khi
  
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  if(slope_type==0)then
     dq=zero
     dbf=zero
     return
  end if

#if NDIM==1
  do n = 1, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              if(slope_type==1.or.slope_type==2)then  ! minmod or average
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
              else
                 write(*,*)'Unknown slope type'
                 stop
              end if
           end do
        end do
     end do     
  end do
#endif

#if NDIM==2              
  if(slope_type==1.or.slope_type==2)then  ! minmod or average
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 2d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    dfll = q(l,i-1,j-1,k,n)-q(l,i,j,k,n)
                    dflm = q(l,i-1,j  ,k,n)-q(l,i,j,k,n)
                    dflr = q(l,i-1,j+1,k,n)-q(l,i,j,k,n)
                    dfml = q(l,i  ,j-1,k,n)-q(l,i,j,k,n)
                    dfmm = q(l,i  ,j  ,k,n)-q(l,i,j,k,n)
                    dfmr = q(l,i  ,j+1,k,n)-q(l,i,j,k,n)
                    dfrl = q(l,i+1,j-1,k,n)-q(l,i,j,k,n)
                    dfrm = q(l,i+1,j  ,k,n)-q(l,i,j,k,n)
                    dfrr = q(l,i+1,j+1,k,n)-q(l,i,j,k,n)
                    
                    vmin = min(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                    vmax = max(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                    
                    dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                    dff  = half*(abs(dfx)+abs(dfy))
                    
                    if(dff>zero)then
                       slop = min(one,min(abs(vmin),abs(vmax))/dff)
                    else
                       slop = one
                    endif
                    
                    dlim = slop
                    
                    dq(l,i,j,k,n,1) = dlim*dfx
                    dq(l,i,j,k,n,2) = dlim*dfy

                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type'
     stop
  endif
  ! 1D transverse TVD slopes for face-centered magnetic fields
  ! Bx along direction Y 
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi+1 ! WARNING HERE
           do l = 1, ngrid
              dlft = slope_type*(bf(l,i,j  ,k,1) - bf(l,i,j-1,k,1))
              drgt = slope_type*(bf(l,i,j+1,k,1) - bf(l,i,j  ,k,1))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dbf(l,i,j,k,1,1) = dsgn*min(dlim,abs(dcen))
           end do
        enddo
     end do
  end do
  ! By along direction X
  do k = klo, khi
     do j = jlo, jhi+1 ! WARNING HERE
        do i = ilo, ihi
           do l = 1, ngrid
              dlft = slope_type*(bf(l,i  ,j,k,2) - bf(l,i-1,j,k,2))
              drgt = slope_type*(bf(l,i+1,j,k,2) - bf(l,i  ,j,k,2))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dbf(l,i,j,k,2,1) = dsgn*min(dlim,abs(dcen))
           end do
        enddo
     end do
  end do
#endif

#if NDIM==3
  if(slope_type==1.or.slope_type==2)then  ! minmod or average
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                    drgt = slope_type*(q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 3d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    dflll = q(l,i-1,j-1,k-1,n)-q(l,i,j,k,n)
                    dflml = q(l,i-1,j  ,k-1,n)-q(l,i,j,k,n)
                    dflrl = q(l,i-1,j+1,k-1,n)-q(l,i,j,k,n)
                    dfmll = q(l,i  ,j-1,k-1,n)-q(l,i,j,k,n)
                    dfmml = q(l,i  ,j  ,k-1,n)-q(l,i,j,k,n)
                    dfmrl = q(l,i  ,j+1,k-1,n)-q(l,i,j,k,n)
                    dfrll = q(l,i+1,j-1,k-1,n)-q(l,i,j,k,n)
                    dfrml = q(l,i+1,j  ,k-1,n)-q(l,i,j,k,n)
                    dfrrl = q(l,i+1,j+1,k-1,n)-q(l,i,j,k,n)

                    dfllm = q(l,i-1,j-1,k  ,n)-q(l,i,j,k,n)
                    dflmm = q(l,i-1,j  ,k  ,n)-q(l,i,j,k,n)
                    dflrm = q(l,i-1,j+1,k  ,n)-q(l,i,j,k,n)
                    dfmlm = q(l,i  ,j-1,k  ,n)-q(l,i,j,k,n)
                    dfmmm = q(l,i  ,j  ,k  ,n)-q(l,i,j,k,n)
                    dfmrm = q(l,i  ,j+1,k  ,n)-q(l,i,j,k,n)
                    dfrlm = q(l,i+1,j-1,k  ,n)-q(l,i,j,k,n)
                    dfrmm = q(l,i+1,j  ,k  ,n)-q(l,i,j,k,n)
                    dfrrm = q(l,i+1,j+1,k  ,n)-q(l,i,j,k,n)

                    dfllr = q(l,i-1,j-1,k+1,n)-q(l,i,j,k,n)
                    dflmr = q(l,i-1,j  ,k+1,n)-q(l,i,j,k,n)
                    dflrr = q(l,i-1,j+1,k+1,n)-q(l,i,j,k,n)
                    dfmlr = q(l,i  ,j-1,k+1,n)-q(l,i,j,k,n)
                    dfmmr = q(l,i  ,j  ,k+1,n)-q(l,i,j,k,n)
                    dfmrr = q(l,i  ,j+1,k+1,n)-q(l,i,j,k,n)
                    dfrlr = q(l,i+1,j-1,k+1,n)-q(l,i,j,k,n)
                    dfrmr = q(l,i+1,j  ,k+1,n)-q(l,i,j,k,n)
                    dfrrr = q(l,i+1,j+1,k+1,n)-q(l,i,j,k,n)

                    vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                         &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                         &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                    vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                         &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                         &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)

                    dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                    dfz  = half*(q(l,i,j,k+1,n)-q(l,i,j,k-1,n))
                    dff  = half*(abs(dfx)+abs(dfy)+abs(dfz))

                    if(dff>zero)then
                       slop = min(one,min(abs(vmin),abs(vmax))/dff)
                    else
                       slop = one
                    endif

                    dlim = slop

                    dq(l,i,j,k,n,1) = dlim*dfx
                    dq(l,i,j,k,n,2) = dlim*dfy
                    dq(l,i,j,k,n,3) = dlim*dfz

                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type'
     stop
  endif     
  ! 2D transverse TVD slopes for face-centered magnetic fields
  ! Bx along direction Y and Z
  xslope_type=MIN(slope_type,2)
  if(slope_type==1.or.slope_type==2.or.slope_type==3)then  ! minmod or average
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi+1 ! WARNING HERE
              ! slopes in first coordinate direction
              do l = 1, ngrid
                 dlft = xslope_type*(bf(l,i,j  ,k,1) - bf(l,i,j-1,k,1))
                 drgt = xslope_type*(bf(l,i,j+1,k,1) - bf(l,i,j  ,k,1))
                 dcen = half*(dlft+drgt)/xslope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,1,1) = dsgn*min(dlim,abs(dcen))
              end do
              ! slopes in second coordinate direction
              do l = 1, ngrid
                 dlft = xslope_type*(bf(l,i,j,k  ,1) - bf(l,i,j,k-1,1))
                 drgt = xslope_type*(bf(l,i,j,k+1,1) - bf(l,i,j,k  ,1))
                 dcen = half*(dlft+drgt)/xslope_type
                 dsgn = sign(one,dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,1,2) = dsgn*min(dlim,abs(dcen))
              end do
           end do
        end do
     end do
  endif
  ! By along direction X and Z
  if(slope_type==1.or.slope_type==2.or.slope_type==3)then  ! minmod or average
     do k = klo, khi
        do j = jlo, jhi+1 ! WARNING HERE
           do i = ilo, ihi
              ! slopes in first coordinate direction
              do l = 1, ngrid
                 dlft = xslope_type*(bf(l,i  ,j,k,2) - bf(l,i-1,j,k,2))
                 drgt = xslope_type*(bf(l,i+1,j,k,2) - bf(l,i  ,j,k,2))
                 dcen = half*(dlft+drgt)/xslope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,2,1) = dsgn*min(dlim,abs(dcen))
              end do
              ! slopes in second coordinate direction
              do l = 1, ngrid
                 dlft = xslope_type*(bf(l,i,j,k  ,2) - bf(l,i,j,k-1,2))
                 drgt = xslope_type*(bf(l,i,j,k+1,2) - bf(l,i,j,k  ,2))
                 dcen = half*(dlft+drgt)/xslope_type
                 dsgn = sign(one,dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,2,2) = dsgn*min(dlim,abs(dcen))
              end do
           end do
        end do
     end do
  endif
  ! Bz along direction X and Y
  if(slope_type==1.or.slope_type==2.or.slope_type==3)then  ! minmod or average
     do k = klo, khi+1 ! WARNING HERE
        do j = jlo, jhi
           do i = ilo, ihi
              ! slopes in first coordinate direction
              do l = 1, ngrid
                 dlft = xslope_type*(bf(l,i  ,j,k,3) - bf(l,i-1,j,k,3))
                 drgt = xslope_type*(bf(l,i+1,j,k,3) - bf(l,i  ,j,k,3))
                 dcen = half*(dlft+drgt)/xslope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,3,1) = dsgn*min(dlim,abs(dcen))
              end do
              ! slopes in second coordinate direction
              do l = 1, ngrid
                 dlft = xslope_type*(bf(l,i,j  ,k,3) - bf(l,i,j-1,k,3))
                 drgt = xslope_type*(bf(l,i,j+1,k,3) - bf(l,i,j  ,k,3))
                 dcen = half*(dlft+drgt)/xslope_type
                 dsgn = sign(one,dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,3,2) = dsgn*min(dlim,abs(dcen))
              end do
           end do
        end do
     end do
  endif
#endif
  
end subroutine uslope
