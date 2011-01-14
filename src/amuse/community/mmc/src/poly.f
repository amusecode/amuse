      subroutine poly
      include 'common.h'

Code to construct a polytrope of index n with Ntot particles inside radius rmax
cThe input parameters are interpreted as follows:
c     imodel = 5
c     W0 is the polytropic index (5 for Plummer)
c     rplum is the tidal radius (in the units of this code, in which the core radius is O(1))
cThe polytropic equation used here is
c     y" + 2y'/r = - y**n, y(0) = 1

c      implicit double precision (a-h,o-z)
      external rkqs,derivsp
      dimension ystart(2)
      PARAMETER (MAXSTP=10000,NVARMAX=2,KMAXX=10000)
      double precision dxsav,xp(KMAXX),
     *yp(NVARMAX,KMAXX),n,y1(kmaxx),y2(kmaxx),xloc(kmaxx),y12(kmaxx),
     &     y22(kmaxx),m(kmaxx),r3(kmaxx),mass,mtot,ke,x0(3,nmax),
     &     xdot0(3,nmax),npoly
      COMMON /path/ kmax,kount,dxsav,xp,yp
      common /parameters/ npoly
c      pi = 4.d0*atan(1.d0)
      print*,'check in poly: pi =',pi
c      read (5,*) n,ntot,rmax,idum
c      print*,'#',n,ntot,rmax
      npoly = w0
      ntot = nt
      rmax = rplum
      print*,'check in poly: making ',ntot,'objects'
      if (ntot.gt.nmax) then
         print*,'nmax too small'
         stop
      endif
      kmax = kmaxx
      dxsav = 0.01d0
      nvar = 2
      x1 = 0.01d0
      x2 = rmax
      h1 = 0.01
      hmin = 0
      ystart(1) = 1.d0 - x1**2/6.d0
      ystart(2) = -x1/3.d0
      eps = 1.d-10
c      print*, '#',ystart,nvar,x1,x2,eps,h1,hmin
      call odeintp(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivsp,rkqs)

cPrepare for spline interpolation of radius**3 as a function of mass, and 
cphi as a function of radius
cCompute energies
      pe = 0.d0
      ke = 0.d0
      r3(1) = 0.d0
      m(1) = 0.d0
      xloc(1) = 0.d0
      y1(1) = 1.d0
      do i = 1,kount
         r3(i+1) = xp(i)**3
         m(i+1) = -xp(i)**2*yp(2,i)
         y1(i+1) = yp(1,i)
         xloc(i+1) = xp(i)
c         print*,'#',xp(i),yp(1,i),yp(2,i),m(i+1),r3(i+1)
         if (i.eq.1) then
            pe = -0.6d0*m(2)**2/xloc(2)
         else
            pe = pe - 0.5d0*(m(i)/xloc(i) + m(i+1)/xloc(i+1))*
     &           (m(i+1) - m(i))
         endif
         ke = ke + 1.5d0*0.5d0*(y1(i) + y1(i+1))*(m(i+1) - m(i))/
     &        (npoly+1.d0)
c         print*,'poly: i,m,r3,y1',i,m(i+1),r3(i+1),y1(i+1)
      enddo
      if (kount.eq.kmaxx) then
         print*,'kount too big; increase kmaxx everywhere'
         stop
      endif
      kount = kount + 1
cInterpolate 
      yp1 = 3.d0
      ypn = 3.d0/yp(1,kount-1)**npoly
      call spline(m,r3,kount,yp1,ypn,y2)
      yp1 = 0.d0
      ypn = yp(2,kount-1)
      call spline(xloc,y1,kount,yp1,ypn,y12)
Create the particles
c      print*,'#',(y2(i),i=1,5)
      mtot = m(kount)
      print*,'poly: kount, mtot, outer radius',kount,mtot,xp(kount-1)
c      stop
      fmax = (npoly-1.5d0)**(npoly-1.5d0)/(npoly-0.5d0)**(npoly-0.5d0)
      do i = 1,ntot
         bodyloc = 1.d0/ntot
         mass = mtot*ran2(idum)
c         print*,'calling splint: mass',mass
         call splint(m,r3,y2,kount,mass,radius3)
         radius = radius3**(1.d0/3.d0)
c         print*,'#',mass,mtot,radius
         sinth = -1.d0 + 2.d0*ran2(idum)
         phi = 2.d0*pi*ran2(idum)
         costh = sqrt(1.d0 - sinth**2)
cNow the velocity
         call splint(xloc,y1,y12,kount,radius,y)         
         print*,'calling splint: radius,y',radius,y

         vmax = sqrt(2.d0*y)
 10      continue
         v = ran2(idum)
         f = v**2*(1-v**2)**(npoly-1.5d0)
         if (f.lt.ran2(idum)*fmax) goto 10
         v = vmax*v
         sinthv = -1.d0 + 2.d0*ran2(idum)
         phiv = 2.d0*pi*ran2(idum)
         costhv = sqrt(1.d0 - sinthv**2)
         print*,bodyloc,radius*costh*cos(phi),radius*costh*sin(phi),
     &        radius*sinth,
     &        v*costhv*cos(phiv),v*costhv*sin(phiv),
     &        v*sinthv
         x0(1,i) = radius*costh*cos(phi)
         x0(2,i) = radius*costh*sin(phi)
         x0(3,i) = radius*sinth
         xdot0(1,i) = v*costhv*cos(phiv)
         xdot0(2,i) = v*costhv*sin(phiv)
         xdot0(3,i) = v*sinthv
      enddo
cPrior to virialisation, convert previous pe and ke to mtot = 1      
      pe = pe/mtot**2
      ke = ke/mtot
      xscale = -pe/0.5d0
      vscale = sqrt(0.25d0/ke)
cNow scale positions and velocities
      do i = 1,ntot
         do k = 1,3
            x0(k,i) = x0(k,i)*xscale
            xdot0(k,i) = xdot0(k,i)*vscale
            x(i,k) = x0(k,i)
            xdot(i,k) = xdot0(k,i)
         enddo
      enddo
cOptional check of the final energies
      pe = 0.d0
      ke = 0.d0
      do i = 1,ntot
         ke = ke + 0.5d0*bodyloc*(xdot0(1,i)**2 + xdot0(2,i)**2 + 
     &        xdot0(3,i)**2)
C          if (i.gt.1) then
C             do j = 1,i-1
C                r2 = 0.d0
C                do k = 1,3
C                   r2 = r2 + (x0(k,i) - x0(k,j))**2
C                enddo
C                pe = pe - bodyloc*bodyloc/sqrt(r2)
C             enddo
C          endif
      enddo
      print*,'N-body system pe (may not be calculated), ke:',pe,ke
CNow set radii needed by MC code
      rtidkg = rplum*xscale
      rtid = rtidkg
c      stop
      end

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
c      REAL x,y,xa(n),y2a(n),ya(n)
      double precision x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
c      REAL a,b,h
      double precision a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.


C       FUNCTION ran2(idum)
C       INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
C c      REAL ran2,AM,EPS,RNMX
C       double precision ran2,AM,EPS,RNMX
C       PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
C      *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
C      *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
C       INTEGER idum2,j,k,iv(NTAB),iy
C       SAVE iv,iy,idum2
C       DATA idum2/123456789/, iv/NTAB*0/, iy/0/
C       if (idum.le.0) then
C         idum=max(-idum,1)
C         idum2=idum
C         do 11 j=NTAB+8,1,-1
C           k=idum/IQ1
C           idum=IA1*(idum-k*IQ1)-k*IR1
C           if (idum.lt.0) idum=idum+IM1
C           if (j.le.NTAB) iv(j)=idum
C 11      continue
C         iy=iv(1)
C       endif
C       k=idum/IQ1
C       idum=IA1*(idum-k*IQ1)-k*IR1
C       if (idum.lt.0) idum=idum+IM1
C       k=idum2/IQ2
C       idum2=IA2*(idum2-k*IQ2)-k*IR2
C       if (idum2.lt.0) idum2=idum2+IM2
C       j=1+iy/NDIV
C       iy=iv(j)-idum2
C       iv(j)=idum
C       if(iy.lt.1)iy=iy+IMM1
C       ran2=min(AM*iy,RNMX)
C       return
C       END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.



      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
c      REAL yp1,ypn,x(n),y(n),y2(n)
      double precision yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=10001)
      INTEGER i,k
c      REAL p,qn,sig,un,u(NMAX)
      double precision p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.


       SUBROUTINE odeintp(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *rkqs)
       INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NVARMAX
c      REAL eps,h1,hmin,x1,x2,ystart(nvar),TINY
       double precision eps,h1,hmin,x1,x2,ystart(nvar),TINY
       EXTERNAL derivs,rkqs
       PARAMETER (MAXSTP=10000,NVARMAX=2,KMAXX=10000,TINY=1.e-30)
       INTEGER i,kmax,kount,nstp
c      REAL dxsav,h,hdid,hnext,x,xsav,dydx(NVARMAX),xp(KMAXX),
c     y(NVARMAX),
       double precision dxsav,h,hdid,hnext,x,xsav,dydx(NVARMAX),
     &     xp(KMAXX),
     &     y(NVARMAX),
     *yp(NVARMAX,KMAXX),yscal(NVARMAX)
       COMMON /path/ kmax,kount,dxsav,xp,yp
       x=x1
       h=sign(h1,x2-x1)
       nok=0
       nbad=0
       kount=0
       do 11 i=1,nvar
         y(i)=ystart(i)
 11    continue
       if (kmax.gt.0) xsav=x-2.*dxsav
       do 16 nstp=1,MAXSTP
         call derivs(x,y,dydx)
         do 12 i=1,nvar
           yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
 12      continue
         if(kmax.gt.0)then
           if(abs(x-xsav).gt.abs(dxsav)) then
             if(kount.lt.kmax-1)then
               kount=kount+1
               xp(kount)=x
               do 13 i=1,nvar
                 yp(i,kount)=y(i)
 13            continue
               xsav=x
             endif
           endif
         endif
         if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
         call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
         if(hdid.eq.h)then
           nok=nok+1
         else
           nbad=nbad+1
         endif
         if((x-x2)*(x2-x1).ge.0.)then
           do 14 i=1,nvar
             ystart(i)=y(i)
 14        continue
           if(kmax.ne.0)then
             kount=kount+1
             xp(kount)=x
             do 15 i=1,nvar
               yp(i,kount)=y(i)
 15          continue
           endif
           return
         endif
         if(abs(hnext).lt.hmin) pause
     *'stepsize smaller than minimum in odeint'
         h=hnext
 16    continue
       pause 'too many steps in odeint'
       return
       END
c  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.

      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER n,NVARMAX
c      REAL eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NVARMAX=2)
CU    USES derivs,rkck
      INTEGER i
C       REAL errmax,h,xnew,yerr(NVARMAX),ytemp(NVARMAX),SAFETY,PGROW,PSHRNK,
       double precision errmax,h,xnew,yerr(NVARMAX),ytemp(NVARMAX),SAFETY,
     &PGROW,PSHRNK,
     *ERRCON

      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        h=SAFETY*h*(errmax**PSHRNK)
        if(h.lt.0.1*h)then
          h=.1*h
        endif
        xnew=x+h
        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.

      subroutine derivsp(x,y,dydx)
      implicit double precision (a-h,o-z)
      dimension y(2),dydx(2)
      double precision npoly
      common /parameters/ npoly
      dydx(1) = y(2)
      dydx(2) = -2*y(2)/x - y(1)**npoly
      return
      end

      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      INTEGER n,NVARMAX
c      REAL h,x,dydx(n),y(n),yerr(n),yout(n)
      double precision h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NVARMAX=2)
CU    USES derivs
      INTEGER i
c      REAL ak2(NVARMAX),ak3(NVARMAX),ak4(NVARMAX),ak5(NVARMAX),ak6(NVARMAX),
      double precision ak2(NVARMAX),ak3(NVARMAX),ak4(NVARMAX),
     &     ak5(NVARMAX),
     &     ak6(NVARMAX),
     *ytemp(NVARMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     *B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     *B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     *B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     *C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     *DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     *DC6=C6-.25)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.
