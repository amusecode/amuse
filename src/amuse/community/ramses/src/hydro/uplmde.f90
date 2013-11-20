!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine tracex(q,dq,c,qm,qp,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer::ngrid
  real(dp)::dx, dt  

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c  

  ! Local variables
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::i, j, k, l, n
  integer ::ir, iu, ip
  real(dp)::dtdx,project_out
  real(dp)::cc, ccc, csq, r, u, p, a
  real(dp)::drx, dux, dpx, dax
  real(dp)::alpham, alphap, alpha0r
  real(dp)::spminus, spplus, spzero
  real(dp)::apright, amright, azrright, azaright
  real(dp)::apleft,  amleft,  azrleft,  azaleft
  
  dtdx = dt/dx
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; ip=3
  project_out=one !zero

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Cell centered values
              cc  = c (l,i,j,k)
              r   = q (l,i,j,k,ir)
              u   = q (l,i,j,k,iu)
              p   = q (l,i,j,k,ip)
              csq = gamma*p/r

              ! TVD slopes in X direction
              drx = dq(l,i,j,k,ir,1)
              dux = dq(l,i,j,k,iu,1)
              dpx = dq(l,i,j,k,ip,1)
              
              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dux) > three*cc)ccc=zero

              ! Characteristic analysis along X direction
              alpham  = half*(dpx/csq - dux*r/cc)
              alphap  = half*(dpx/csq + dux*r/cc)
              alpha0r = drx - dpx/csq

              ! Right state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)>zero)spplus =-project_out
              if((u-ccc)>zero)spminus=-project_out
              if( u     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap
              amright  = half*(-one-spminus)*alpham 
              azrright = half*(-one-spzero )*alpha0r

              qp(l,i,j,k,ir,1) = r + (apright+amright+azrright)
              qp(l,i,j,k,iu,1) = u + (apright-amright         )*cc/r
              qp(l,i,j,k,ip,1) = p + (apright+amright         )*csq
              qp(l,i,j,k,ir,1) = max(smallr,qp(l,i,j,k,ir,1))

              ! Left state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)<=zero)spplus =+project_out
              if((u-ccc)<=zero)spminus=+project_out
              if( u     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap
              amleft   = half*(+one-spminus)*alpham
              azrleft  = half*(+one-spzero )*alpha0r

              qm(l,i,j,k,ir,1) = r + (apleft+amleft+azrleft)
              qm(l,i,j,k,iu,1) = u + (apleft-amleft        )*cc/r
              qm(l,i,j,k,ip,1) = p + (apleft+amleft        )*csq
              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))

           end do
        end do
     end do
  end do

#if NVAR > NDIM + 2
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a   =  q(l,i,j,k,n)    ! Cell centered values
                 u   =  q(l,i,j,k,iu)
                 dax = dq(l,i,j,k,n,1)  ! TVD slopes
                 
                 ! Right state
                 spzero=(u    )*dtdx
                 if(u>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*dax
                 qp(l,i,j,k,n,1) = a + azaright

                 ! Left state
                 spzero=(u    )*dtdx
                 if(u<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero)*dax
                 qm(l,i,j,k,n,1) = a + azaleft
              end do
           end do
        end do
     end do
  end do
#endif

end subroutine tracex
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>1
subroutine tracexy(q,dq,c,qm,qp,dx,dy,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c  

  ! declare local variables
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::i, j, k, l, n
  integer ::ir, iu, iv, ip
  real(dp)::dtdx,dtdy,project_out
  real(dp)::cc, ccc, csq, r, u, v, p, a
  real(dp)::drx, dux, dvx, dpx, dax
  real(dp)::dry, duy, dvy, dpy, day
  real(dp)::alpham, alphap, alpha0r, alpha0u, alpha0v
  real(dp)::spminus, spplus, spzero
  real(dp)::apright, amright, azrright, azuright, azvright, azaright
  real(dp)::apleft,  amleft,  azrleft,  azuleft,  azvleft,  azaleft
  real(dp)::srx,sux,svx,spx,sax
  real(dp)::sry,suy,svy,spy,say
    
  dtdx = dt/dx; dtdy = dt/dy
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; ip=4
  project_out=one !zero

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! cell centered values
              cc  = c  (l,i,j,k)
              r   = q  (l,i,j,k,ir)
              u   = q  (l,i,j,k,iu)
              v   = q  (l,i,j,k,iv)
              p   = q  (l,i,j,k,ip)
              csq = gamma*p/r

              ! TVD slopes in X and Y directions
              drx = dq(l,i,j,k,ir,1)
              dux = dq(l,i,j,k,iu,1)
              dvx = dq(l,i,j,k,iv,1)
              dpx = dq(l,i,j,k,ip,1)
              
              dry = dq(l,i,j,k,ir,2)
              duy = dq(l,i,j,k,iu,2)
              dvy = dq(l,i,j,k,iv,2)
              dpy = dq(l,i,j,k,ip,2)
              
              ! Transverse derivatives
              srx = half*dtdy*(-v*dry - (dvy)*r      )
              sux = half*dtdy*(-v*duy                )
              svx = half*dtdy*(-v*dvy - (dpy)/r      )
              spx = half*dtdy*(-v*dpy - (dvy)*gamma*p)

              sry = half*dtdx*(-u*drx - (dux)*r      )
              suy = half*dtdx*(-u*dux - (dpx)/r      )
              svy = half*dtdx*(-u*dvx                )
              spy = half*dtdx*(-u*dpx - (dux)*gamma*p)

              ! Characteristic analysis along X direction
              alpham  = half*(dpx/csq - dux*r/cc)
              alphap  = half*(dpx/csq + dux*r/cc)
              alpha0r = drx - dpx/csq
              alpha0v = dvx

              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dux) > three*cc)ccc=zero

              ! Right state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)>zero)spplus =-project_out
              if((u-ccc)>zero)spminus=-project_out
              if( u     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap 
              amright  = half*(-one-spminus)*alpham 
              azrright = half*(-one-spzero )*alpha0r
              azvright = half*(-one-spzero )*alpha0v

              qp(l,i,j,k,ir,1) = r + (apright+amright+azrright)     +srx
              qp(l,i,j,k,iu,1) = u + (apright-amright         )*cc/r+sux
              qp(l,i,j,k,ip,1) = p + (apright+amright         )*csq +spx
              qp(l,i,j,k,iv,1) = v + (                azvright)     +svx
              qp(l,i,j,k,ir,1) = max(smallr,qp(l,i,j,k,ir,1))

              ! Left state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)<=zero)spplus =+project_out
              if((u-ccc)<=zero)spminus=+project_out
              if( u     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap 
              amleft   = half*(+one-spminus)*alpham 
              azrleft  = half*(+one-spzero )*alpha0r
              azvleft  = half*(+one-spzero )*alpha0v
              
              qm(l,i,j,k,ir,1) = r + (apleft+amleft+azrleft)     +srx
              qm(l,i,j,k,iu,1) = u + (apleft-amleft        )*cc/r+sux
              qm(l,i,j,k,ip,1) = p + (apleft+amleft        )*csq +spx
              qm(l,i,j,k,iv,1) = v + (              azvleft)     +svx
              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))

              ! Characteristic analysis along Y direction
              alpham  = half*(dpy/csq - dvy*r/cc)
              alphap  = half*(dpy/csq + dvy*r/cc)
              alpha0r = dry - dpy/csq
              alpha0u = duy

              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dvy) > three*cc)ccc=zero

              ! Top state
              spminus = (v-ccc)*dtdy
              spplus  = (v+ccc)*dtdy
              spzero  = (v    )*dtdy
              if((v+ccc)>zero)spplus =-project_out
              if((v-ccc)>zero)spminus=-project_out
              if( v     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap
              amright  = half*(-one-spminus)*alpham 
              azrright = half*(-one-spzero )*alpha0r
              azuright = half*(-one-spzero )*alpha0u

              qp(l,i,j,k,ir,2) = r + (apright+amright+azrright)     +sry
              qp(l,i,j,k,iv,2) = v + (apright-amright         )*cc/r+svy
              qp(l,i,j,k,ip,2) = p + (apright+amright         )*csq +spy
              qp(l,i,j,k,iu,2) = u + (                azuright)     +suy
              qp(l,i,j,k,ir,2) = max(smallr,qp(l,i,j,k,ir,2))

              ! Bottom state
              spminus = (v-ccc)*dtdy
              spplus  = (v+ccc)*dtdy
              spzero  = (v    )*dtdy
              if((v+ccc)<=zero)spplus =+project_out
              if((v-ccc)<=zero)spminus=+project_out
              if( v     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap
              amleft   = half*(+one-spminus)*alpham
              azrleft  = half*(+one-spzero )*alpha0r
              azuleft  = half*(+one-spzero )*alpha0u

              qm(l,i,j,k,ir,2) = r + (apleft+amleft+azrleft)     +sry
              qm(l,i,j,k,iv,2) = v + (apleft-amleft        )*cc/r+svy
              qm(l,i,j,k,ip,2) = p + (apleft+amleft        )*csq +spy
              qm(l,i,j,k,iu,2) = u + (              azuleft)     +suy
              qm(l,i,j,k,ir,2) = max(smallr, qm(l,i,j,k,ir,2))

           end do
        end do
     end do
  end do

#if NVAR > NDIM + 2
  ! Passive scalars
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a = q(l,i,j,k,n)     ! Cell centered values
                 u = q(l,i,j,k,iu)
                 v = q(l,i,j,k,iv)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 sax = half*dtdy*(-v*day) ! Transverse
                 say = half*dtdx*(-u*dax) ! derivatives

                 ! Right state
                 spzero=(u    )*dtdx
                 if(u>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*dax
                 qp(l,i,j,k,n,1) = a + azaright + sax
                 
                 ! Left state
                 spzero=(u    )*dtdx
                 if(u<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*dax
                 qm(l,i,j,k,n,1) = a + azaleft + sax
                 
                 ! Top state
                 spzero=(v    )*dtdy
                 if(v>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*day
                 qp(l,i,j,k,n,2) = a + azaright + say
                 
                 ! Bottom state
                 spzero=(v    )*dtdy
                 if(v<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*day
                 qm(l,i,j,k,n,2) = a + azaleft + say
                 
              end do
           end do
        end do
     end do
  end do
#endif

end subroutine tracexy
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>2
subroutine tracexyz(q,dq,c,qm,qp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx,dy,dz, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c  

  ! declare local variables
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::i, j, k, l, n
  integer ::ir, iu, iv, iw, ip
  real(dp)::dtdx,dtdy,dtdz,project_out
  real(dp)::cc, ccc, csq, r, u, v, w, p, a
  real(dp)::drx, dux, dvx, dwx, dpx, dax
  real(dp)::dry, duy, dvy, dwy, dpy, day
  real(dp)::drz, duz, dvz, dwz, dpz, daz
  real(dp)::alpham, alphap, alpha0r, alpha0u, alpha0v, alpha0w
  real(dp)::spminus, spplus, spzero
  real(dp)::apright, amright, azrright, azuright, azvright, azwright, azaright
  real(dp)::apleft,  amleft,  azrleft,  azuleft,  azvleft,  azwleft,  azaleft
  real(dp)::srx,sux,svx,swx,spx,sax
  real(dp)::sry,suy,svy,swy,spy,say
  real(dp)::srz,suz,svz,swz,spz,saz
    
  dtdx = dt/dx; dtdy = dt/dy; dtdz = dt/dz
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5
  project_out=one !zero
  
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
  
              ! Cell centered values
              cc  = c  (l,i,j,k)
              r   = q  (l,i,j,k,ir)
              u   = q  (l,i,j,k,iu)
              v   = q  (l,i,j,k,iv)
              w   = q  (l,i,j,k,iw)
              p   = q  (l,i,j,k,ip)
              csq = gamma*p/r

              ! TVD slopes in all 3 directions
              drx = dq(l,i,j,k,ir,1)
              dux = dq(l,i,j,k,iu,1)
              dvx = dq(l,i,j,k,iv,1)
              dwx = dq(l,i,j,k,iw,1)
              dpx = dq(l,i,j,k,ip,1)
              
              dry = dq(l,i,j,k,ir,2)
              duy = dq(l,i,j,k,iu,2)
              dvy = dq(l,i,j,k,iv,2)
              dwy = dq(l,i,j,k,iw,2)
              dpy = dq(l,i,j,k,ip,2)
              
              drz = dq(l,i,j,k,ir,3)
              duz = dq(l,i,j,k,iu,3)
              dvz = dq(l,i,j,k,iv,3)
              dwz = dq(l,i,j,k,iw,3)
              dpz = dq(l,i,j,k,ip,3)

              ! Transverse derivatives
              srx = half*dtdx*(-v*dry-w*drz - (dvy+dwz)*r      )
              spx = half*dtdx*(-v*dpy-w*dpz - (dvy+dwz)*gamma*p)
              sux = half*dtdx*(-v*duy-w*duz                    )
              svx = half*dtdx*(-v*dvy-w*dvz - (dpy)/r          )
              swx = half*dtdx*(-v*dwy-w*dwz - (dpz)/r          )

              sry = half*dtdx*(-u*drx-w*drz - (dux+dwz)*r      )
              spy = half*dtdx*(-u*dpx-w*dpz - (dux+dwz)*gamma*p)
              suy = half*dtdx*(-u*dux-w*duz - (dpx)/r          )
              svy = half*dtdx*(-u*dvx-w*dvz                    )
              swy = half*dtdx*(-u*dwx-w*dwz - (dpz)/r          )

              srz = half*dtdx*(-v*dry-u*drx - (dvy+dux)*r      )
              spz = half*dtdx*(-v*dpy-u*dpx - (dvy+dux)*gamma*p)
              suz = half*dtdx*(-v*duy-u*dux - (dpx)/r          )
              svz = half*dtdx*(-v*dvy-u*dvx - (dpy)/r          )
              swz = half*dtdx*(-v*dwy-u*dwx                    )

              ! Characteristic analysis along X direction
              alpham  = half*(dpx/csq - dux*r/cc)
              alphap  = half*(dpx/csq + dux*r/cc)
              alpha0r = drx - dpx/csq
              alpha0v = dvx
              alpha0w = dwx

              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dux) > three*cc)ccc=zero

              ! Right state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)>zero)spplus =-project_out
              if((u-ccc)>zero)spminus=-project_out
              if( u     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap
              amright  = half*(-one-spminus)*alpham
              azrright = half*(-one-spzero )*alpha0r
              azvright = half*(-one-spzero )*alpha0v
              azwright = half*(-one-spzero )*alpha0w

              qp(l,i,j,k,ir,1) = r + (apright+amright+azrright)     +srx
              qp(l,i,j,k,iu,1) = u + (apright-amright         )*cc/r+sux
              qp(l,i,j,k,ip,1) = p + (apright+amright         )*csq +spx
              qp(l,i,j,k,iv,1) = v + (                azvright)     +svx
              qp(l,i,j,k,iw,1) = w + (                azwright)     +swx
              qp(l,i,j,k,ir,1) = max(smallr,qp(l,i,j,k,ir,1))

              ! Left state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)<=zero)spplus =+project_out
              if((u-ccc)<=zero)spminus=+project_out
              if( u     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap
              amleft   = half*(+one-spminus)*alpham
              azrleft  = half*(+one-spzero )*alpha0r
              azvleft  = half*(+one-spzero )*alpha0v
              azwleft  = half*(+one-spzero )*alpha0w

              qm(l,i,j,k,ir,1) = r + (apleft+amleft+azrleft)     +srx
              qm(l,i,j,k,iu,1) = u + (apleft-amleft        )*cc/r+sux
              qm(l,i,j,k,ip,1) = p + (apleft+amleft        )*csq +spx
              qm(l,i,j,k,iv,1) = v + (              azvleft)     +svx
              qm(l,i,j,k,iw,1) = w + (              azwleft)     +swx
              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))

              ! Characteristic analysis along Y direction
              alpham  = half*(dpy/csq - dvy*r/cc)
              alphap  = half*(dpy/csq + dvy*r/cc)
              alpha0r = dry - dpy/csq
              alpha0u = duy
              alpha0w = dwy

              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dvy) > three*cc)ccc=zero

              ! Top state
              spminus = (v-ccc)*dtdy
              spplus  = (v+ccc)*dtdy
              spzero  = (v    )*dtdy
              if((v+ccc)>zero)spplus =-project_out
              if((v-ccc)>zero)spminus=-project_out
              if( v     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap 
              amright  = half*(-one-spminus)*alpham 
              azrright = half*(-one-spzero )*alpha0r
              azuright = half*(-one-spzero )*alpha0u
              azwright = half*(-one-spzero )*alpha0w

              qp(l,i,j,k,ir,2) = r + (apright+amright+azrright)     +sry
              qp(l,i,j,k,iv,2) = v + (apright-amright         )*cc/r+svy
              qp(l,i,j,k,ip,2) = p + (apright+amright         )*csq +spy
              qp(l,i,j,k,iu,2) = u + (                azuright)     +suy
              qp(l,i,j,k,iw,2) = w + (                azwright)     +swy
              qp(l,i,j,k,ir,2) = max(smallr,qp(l,i,j,k,ir,2))

              ! Bottom state
              spminus = (v-ccc)*dtdy
              spplus  = (v+ccc)*dtdy
              spzero  = (v    )*dtdy
              if((v+ccc)<=zero)spplus =+project_out
              if((v-ccc)<=zero)spminus=+project_out
              if( v     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap 
              amleft   = half*(+one-spminus)*alpham 
              azrleft  = half*(+one-spzero )*alpha0r
              azuleft  = half*(+one-spzero )*alpha0u
              azwleft  = half*(+one-spzero )*alpha0w

              qm(l,i,j,k,ir,2) = r + (apleft+amleft+azrleft)     +sry
              qm(l,i,j,k,iv,2) = v + (apleft-amleft        )*cc/r+svy
              qm(l,i,j,k,ip,2) = p + (apleft+amleft        )*csq +spy
              qm(l,i,j,k,iu,2) = u + (              azuleft)     +suy
              qm(l,i,j,k,iw,2) = w + (              azwleft)     +swy
              qm(l,i,j,k,ir,2) = max(smallr, qm(l,i,j,k,ir,2))

              ! Characteristic analysis along Z direction
              alpham  = half*(dpz/csq - dwz*r/cc)
              alphap  = half*(dpz/csq + dwz*r/cc)
              alpha0r = drz - dpz/csq
              alpha0u = duz
              alpha0v = dvz

              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dwz) > three*cc)ccc=zero

              ! Front state
              spminus = (w-ccc)*dtdz
              spplus  = (w+ccc)*dtdz
              spzero  = (w    )*dtdz
              if((w+ccc)>zero)spplus =-project_out
              if((w-ccc)>zero)spminus=-project_out
              if( w     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap 
              amright  = half*(-one-spminus)*alpham 
              azrright = half*(-one-spzero )*alpha0r
              azuright = half*(-one-spzero )*alpha0u
              azvright = half*(-one-spzero )*alpha0v

              qp(l,i,j,k,ir,3) = r + (apright+amright+azrright)     +srz
              qp(l,i,j,k,iw,3) = w + (apright-amright         )*cc/r+swz
              qp(l,i,j,k,ip,3) = p + (apright+amright         )*csq +spz
              qp(l,i,j,k,iu,3) = u + (                azuright)     +suz
              qp(l,i,j,k,iv,3) = v + (                azvright)     +svz
              qp(l,i,j,k,ir,3) = max(smallr,qp(l,i,j,k,ir,3))

              ! Back state
              spminus = (w-ccc)*dtdz
              spplus  = (w+ccc)*dtdz
              spzero  = (w    )*dtdz
              if((w+ccc)<=zero)spplus =+project_out
              if((w-ccc)<=zero)spminus=+project_out
              if( w     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap 
              amleft   = half*(+one-spminus)*alpham 
              azrleft  = half*(+one-spzero )*alpha0r
              azuleft  = half*(+one-spzero )*alpha0u
              azvleft  = half*(+one-spzero )*alpha0v

              qm(l,i,j,k,ir,3) = r + (apleft+amleft+azrleft)     +srz
              qm(l,i,j,k,iw,3) = w + (apleft-amleft        )*cc/r+swz
              qm(l,i,j,k,ip,3) = p + (apleft+amleft        )*csq +spz
              qm(l,i,j,k,iu,3) = u + (              azuleft)     +suz
              qm(l,i,j,k,iv,3) = v + (              azvleft)     +svz
              qm(l,i,j,k,ir,3) = max(smallr, qm(l,i,j,k,ir,3))
           end do

        end do
     end do
  end do

#if NVAR > NDIM + 2
  ! Passive scalars
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a   =  q(l,i,j,k,n)    ! Cell centered values
                 u   =  q(l,i,j,k,iu)
                 v   =  q(l,i,j,k,iv)
                 w   =  q(l,i,j,k,iw)
                 dax = dq(l,i,j,k,n,1)  ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 daz = dq(l,i,j,k,n,3)
                 sax = half*dtdx*(-v*day-w*daz) ! Transverse
                 say = half*dtdx*(-u*dax-w*daz) ! derivatives
                 saz = half*dtdx*(-v*day-u*dax) ! 

                 
                 ! Right state
                 spzero = (u    )*dtdx
                 if(u>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*dax
                 qp(l,i,j,k,n,1) = a + azaright + sax
                 
                 ! Left state
                 spzero = (u    )*dtdx
                 if(u<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*dax
                 qm(l,i,j,k,n,1) = a + azaleft + sax

                 ! Top state
                 spzero = (v    )*dtdy
                 if(v>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*day
                 qp(l,i,j,k,n,2) = a + azaright + say
                 
                 ! Bottom state
                 spzero = (v    )*dtdy
                 if(v<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*day
                 qm(l,i,j,k,n,2) = a + azaleft + say

                 ! Front state
                 spzero = (w    )*dtdy
                 if(w>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*daz
                 qp(l,i,j,k,n,3) = a + azaright + saz
                 
                 ! Back state
                 spzero = (w    )*dtdy
                 if(w<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*daz
                 qm(l,i,j,k,n,3) = a + azaleft + saz
              end do
           end do
        end do
     end do
  end do
#endif

end subroutine tracexyz
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdivu(q,div,dx,dy,dz,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::div

  integer::i, j, k, l
  real(dp)::factorx, factory, factorz
  real(dp),dimension(1:nvector)::ux, vy, wz

  factorx=half**(ndim-1)/dx
  factory=half**(ndim-1)/dy
  factorz=half**(ndim-1)/dz

  do k = kf1, kf2
     do j = jf1, jf2
        do i = if1, if2

           ux = zero; vy=zero; wz=zero

           if(ndim>0)then
              do l=1, ngrid
                 ux(l)=ux(l)+factorx*(q(l,i  ,j  ,k  ,2) - q(l,i-1,j  ,k  ,2))
              end do
           end if

#if NDIM>1           
           if(ndim>1)then
              do l=1, ngrid
                 ux(l)=ux(l)+factorx*(q(l,i  ,j-1,k  ,2) - q(l,i-1,j-1,k  ,2))
                 vy(l)=vy(l)+factory*(q(l,i  ,j  ,k  ,3) - q(l,i  ,j-1,k  ,3)+&
                      &               q(l,i-1,j  ,k  ,3) - q(l,i-1,j-1,k  ,3))
              end do
           end if
#endif
#if NDIM>2
           if(ndim>2)then
              do l=1, ngrid
                 ux(l)=ux(l)+factorx*(q(l,i  ,j  ,k-1,2) - q(l,i-1,j  ,k-1,2)+&
                      &               q(l,i  ,j-1,k-1,2) - q(l,i-1,j-1,k-1,2))
                 vy(l)=vy(l)+factory*(q(l,i  ,j  ,k-1,3) - q(l,i  ,j-1,k-1,3)+&
                      &               q(l,i-1,j  ,k-1,3) - q(l,i-1,j-1,k-1,3))
                 wz(l)=wz(l)+factorz*(q(l,i  ,j  ,k  ,4) - q(l,i  ,j  ,k-1,4)+&
                      &               q(l,i  ,j-1,k  ,4) - q(l,i  ,j-1,k-1,4)+&
                      &               q(l,i-1,j  ,k  ,4) - q(l,i-1,j  ,k-1,4)+&
                      &               q(l,i-1,j-1,k  ,4) - q(l,i-1,j-1,k-1,4))
              end do
           end if
#endif
           do l=1,ngrid
              div(l,i,j,k) = ux(l) + vy(l) + wz(l)
           end do

        end do
     end do
  end do

end subroutine cmpdivu
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine consup(uin,flux,div,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin 
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::div 

  integer:: i, j, k, l, n
  real(dp)::factor
  real(dp),dimension(1:nvector),save:: div1

  factor=half**(ndim-1)

  ! Add diffusive flux where flow is compressing
  do n = 1, nvar

     do k = kf1, MAX(kf1,ku2-2)
        do j = jf1, MAX(jf1, ju2-2) 
           do i = if1, if2
              div1 = zero
              do l = 1, ngrid
                 div1(l) = factor*div(l,i,j,k)
              end do
#if NDIM>1
              do l = 1, ngrid
                 div1(l)=div1(l)+factor*div(l,i,j+1,k)
              end do
#endif
#if NDIM>2
              do l = 1, ngrid
                 div1(l)=div1(l)+factor*(div(l,i,j,k+1)+div(l,i,j+1,k+1))
              end do
#endif
              do l = 1, ngrid
                 div1(l) = difmag*min(zero,div1(l))
              end do
              do l = 1, ngrid
                 flux(l,i,j,k,n,1) = flux(l,i,j,k,n,1) + &
                      &  dt*div1(l)*(uin(l,i,j,k,n) - uin(l,i-1,j,k,n))
              end do

           end do
        end do
     end do

#if NDIM>1
     do k = kf1, MAX(kf1,ku2-2)
        do j = jf1, jf2
           do i = iu1+2, iu2-2
              div1 = zero
              do l = 1, ngrid
                 div1(l)=div1(l)+factor*(div(l,i,j,k ) + div(l,i+1,j,k))
              end do
#if NDIM>2
              do l = 1, ngrid
                 div1(l)=div1(l)+factor*(div(l,i,j,k+1) + div(l,i+1,j,k+1))
              end do
#endif
              do l = 1, ngrid
                 div1(l) = difmag*min(zero,div1(l))
              end do
              do l = 1, ngrid
                 flux(l,i,j,k,n,2) = flux(l,i,j,k,n,2) + &
                      &  dt*div1(l)*(uin(l,i,j,k,n) - uin(l,i,j-1,k,n))
              end do
           end do
        end do
     end do
#endif

#if NDIM>2
     do k = kf1, kf2
        do j = ju1+2, ju2-2 
           do i = iu1+2, iu2-2 
              do l = 1, ngrid
                 div1(l)=factor*(div(l,i,j  ,k) + div(l,i+1,j  ,k) &
                      &        + div(l,i,j+1,k) + div(l,i+1,j+1,k))
              end do
              do l = 1, ngrid
                 div1(l) = difmag*min(zero,div1(l))
              end do
              do l = 1, ngrid
                 flux(l,i,j,k,n,3) = flux(l,i,j,k,n,3) + &
                      &  dt*div1(l)*(uin(l,i,j,k,n) - uin(l,i,j,k-1,n))
              end do
           end do
        end do
     end do
#endif

  end do

end subroutine consup
