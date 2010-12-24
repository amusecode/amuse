*
*
c      subroutine newpos(k,nup,rnew1,rnew2,iruntemp)
      subroutine newpos(k,nup)
*
*
*       determine new positions of two interacting objects after encounter
*       ------------------------------------------------------------------
*       M. HENON 1971, Astrophysics and Space Science Vol. 14, 151
*       ----------------------------------------------------------
*
*
      include 'common.h'
*
*
      real*8 e1,e2,a1,a2,q1,dz,b1,c1,gmin1,gmax1,gmin2,gmax2,
     &       rmin1,rmin2,rmax1,rmax2,rmaxz,s,u1,u2,v1,drds,q0,
     &       rtot,rles,xcut,ycut,xrr,rnn1,rnn2,ytot,qmax,rnew1,
     &       rnew2
c     &       rtot,rles,xcut,ycut,xrr,vrp,vtp,rnn1,rnn2
*
      real*4 ran2
*
      integer k,im1,im2,n,i,ipoin1,ipoin2,irm1,irm2,ipp,nup,m,imi1,
     &        imi2,ipo,kmin1,kmax1,kmin2,kmax2,imk1,imk2,maxk,
     &        mmin,mmax
*
c      common /timep/ vrp(nmax),vtp(nmax)
*
*
      n = nt
      i = 0
      im1 = iname(k)
      im2 = iname(k+1)
      gmin1 = 0.0d0
      gmin2 = 0.0d0
      gmax1 = 0.0d0
      gmax2 = 0.0d0
      u1 = 0.d0
      u2 = 0.d0
      v1 = 0.d0
      q0 = 0.d0
      ipoin1 = 0
      ipoin2 = 0
      irm1 = 0
      irm2 = 0
      imi1 = 0
      imi2 = 0
      imk1 = 0
      imk2 = 0
c      xmin(im1) = 0.0d0
c      xmin(im2) = 0.0d0
c      xmax(im1) = 0.0d0
c      xmax(im2) = 0.0d0
      xmaxz(im1) = 0.0d0
      xmaxz(im2) = 0.0d0
      xgmin(im1) = 0.0d0
      xgmin(im2) = 0.0d0
      xgmax(im1) = 0.0d0
      xgmax(im2) = 0.0d0
      rnew1 = r(k)
      rnew2 = r(k+1)
      vrp(im1) = vr(im1)
      vtp(im1) = vt(im1)
      vrp(im2) = vr(im2)
      vtp(im2) = vt(im2)
*
*       find maximal distance for new position of objects in super-zone
*
      ipp = nzst(nup)
      maxk = ipp + 1
      
*
      if(ipp.eq.n) then
        rmaxz = 1.d12
        maxk = n
      else
        rmaxz = 0.01d0*r(ipp) + 0.99d0*r(ipp+1)
      endif
*     
*       determine total energy and angular momentum for interacting objects
*
      e1 = u(k) + 0.5d0*(vr(im1)**2 + vt(im1)**2)
      e2 = u(k+1) + 0.5d0*(vr(im2)**2 + vt(im2)**2)
      a1 = r(k)*vt(im1)
      a2 = r(k+1)*vt(im2)
*
*       check for escapers  -   energy > 0
*       do not find a new possition for massles particles
*
      if(e1.gt.0.0d0.or.body(im1).eq.0.d0) then

*
c          xmin(im1) = rnew1
c          xmax(im1) = rnew1
          xmaxz(im1) = rmaxz
          ipoin1 = 1
*
       endif
*
*       do not find a new possition for massles particles
*
      if(e2.gt.0.0d0.or.body(im2).eq.0.d0) then
*
c          xmin(im2) = rnew2
c          xmax(im2) = rnew2
          xmaxz(im2) = rmaxz
          ipoin2 = 1
*
      endif
*
*
*       find rmin and rmax for interacting objects
*
*       first star
*
      if(ipoin1.eq.1) go to 30
*
cThis call attempts to find m for a simple case (the majority).  Otherwise it returns
czero, and then the standard procedure is followed.
      m = mmin(k,e1,a1,n)
      if (m.gt.0) goto 13
      i = -1
   10 i = i + 1
      m = k - i
*
      if(m.eq.0) then
        b1 = 0.0d0
        c1 = u(1)
        go to 15
      endif
*
      q1 = 2.0d0*(e1 - u(m)) - a1**2/r(m)**2
      if(q1.ge.0.0d0) go to 10
*
      if(m.eq.k) then
        imk1 = 1
        rmin1 = 0.9999d0*r(k)
        write(6,*) '< r(k)   k,im1,r(k),u(k),e1,vr,vt,body = ',
     &  k,im1,r(k),u(k),e1,vr(im1),vt(im1),body(im1)
	call flush(6)
        kmin1 = 1
        go to 17
      endif
*
      if(m.eq.n) then
        b1 = u(m)*r(m)
        c1 = 0.0d0
        go to 15
      endif
*
 13   continue
      dz = 1.0d0/r(m+1) - 1.0d0/r(m)
      b1 = (u(m+1) - u(m))/dz
      c1 = (u(m)/r(m+1) - u(m+1)/r(m))/dz
   15 rmin1 = (-b1 + sqrt(b1*b1 - 2.0d0*a1*a1*(c1 - e1)))/a1/a1
      rmin1 = 1.0d0/rmin1
*
      if(rmin1.gt.r(k)) then 
*
        xrr = 1.d0 - rmin1/r(k)
        xrr = abs(xrr)
*
        if(xrr.gt.1.0d-6) then
          write(6,*) '    rmin1  >  r(k)  iseed = ',iseed
          write(6,*) 'm,k,rmin1,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1 =',
     &                m,k,rmin1,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1
          call flush(6)
*
          rmin1 = 0.9999d0*r(k)
          imi1 = 1
        else
          rmin1 = r(k)
          imi1 = 1
        endif
*
      endif
*
      gmin1 = 2.0d0*a1*a1/rmin1**3 + 2.0d0*b1/rmin1**2
c      print*,gmin1,rmin1
*
      if(gmin1.lt.0.0d0) then
*
        write(6,*) '    gmin1  <   0.0   iseed = ',iseed
        write(6,*) 'm,k,rm1,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1,b1,gm1=',
     &            m,k,rmin1,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1,b1,gmin1

        call flush(6)
*
        gmin1 = 1.d10
      endif
*
      if(m.le.1) then
        kmin1 = 1
      else
        kmin1 = m - 1
      endif
*
      if((m.eq.n).or.(imi1.eq.1)) kmin1 = 1
*
 17   continue
*
      m = mmax(k,e1,a1,n)
      if (m.gt.0) goto 23

      i = -1
   20 i = i + 1
      m = k + i
*
      if(m.gt.n) then
        b1 = u(m-1)*r(m-1)
        c1 = 0.0d0
        go to 25
      endif
*
      q1 = 2.0d0*(e1 - u(m)) - a1**2/r(m)**2
c      print*,i,m,q1,e1,u(m),a1,r(m)
      if(q1.ge.0.0d0) go to 20
*
      if(m.eq.k) then
        if(imk1.eq.1) then
          rnew1 = r(k)
          imk1 = 2
          write(6,*) '-- O --  k,im1,r(k),u(k),e1= ',k,im1,r(k),u(k),e1
	  call flush(6)
          rmin1 = rnew1
          rmax1 = rnew1
          go to 30
        endif
*
        imk1 = 1
        rmax1 = 1.0001d0*r(k)
        write(6,*) '> r(k)    k,im1,r(k),u(k),e1 = ',k,im1,r(k),u(k),e1
	call flush(6)
        kmax1 = maxk
        go to 30
      endif
*
      if(m.eq.1) then
*
        if(imi1.eq.1) then
*
        write(6,*) ' m = 1   rnew1 = r(k)  iseed = ',iseed
        write(6,*) 'm,k,r(k),r(m),u(m),e1,a1 =',      
     &              m,k,r(k),r(m),u(m),e1,a1      
        call flush(6)
*
          imi1 = 2
          rnew1 = r(k)
          go to 30
        endif
*
        write(6,*) ' m = 1   rmax1 = 1.0001*r(k)   iseed = ',iseed
        write(6,*) 'm,k,r(k),r(m),u(m),e1,a1 =',      
     &              m,k,r(k),r(m),u(m),e1,a1      
        call flush(6)
*
        rmax1 = 1.0001d0*r(k)
        gmax1 = -1.d10
        go to 27
      endif
*
 23   continue
      dz = 1.0d0/r(m) - 1.0d0/r(m-1)
      b1 = (u(m) - u(m-1))/dz
      c1 = (u(m-1)/r(m) - u(m)/r(m-1))/dz
   25 continue
c      print*,c1,u(m-1),r(m),u(m),r(m-1),dz,m,k
      rmax1 = 0.5d0*(-b1 + sqrt(b1*b1 - 2.0d0*a1*a1*(c1 - e1)))
     &        /(c1 - e1)
*
      if(rmax1.lt.r(k)) then 
*
        if(imi1.eq.1) then
          write(6,*) 'circular orbit     iseed,k = ',iseed,k
	  call flush(6)
          imi1 = 2
          rnew1 = r(k)
          rmin1 = rnew1
          rmax1 = rnew1
          go to 30
        endif
*        
        xrr = 1.d0 - rmax1/r(k)
        xrr = abs(xrr)
*
        if(xrr.gt.1.0d-6) then
          write(6,*) ' rmax1 < r(k)  iseed =',iseed
          write(6,*) 'm,k,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1 =',
     &                m,k,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1          
          call flush(6)
*          
          rmax1 = 1.0001d0*r(k)
          imi1 = 3
        else
          rmax1 = r(k)
          imi1 = 3
        endif
*
      endif
*
      gmax1 = 2.0d0*a1*a1/rmax1**3 + 2.0d0*b1/rmax1**2
*
      if(gmax1.gt.0.0d0) then
*
        write(6,*) ' gmax1  >  0.0   iseed,gmax1 = ',iseed,gmax1
        write(6,*) 'm,k,rmax1,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1 =',
     &              m,k,rmax1,r(k),r(m),r(m-1),u(m),u(m-1),e1,a1          
        call flush(6)
*          
        gmax1 = -1.d10
      endif
*
  27  if(rmax1.lt.rmin1) then
*
        write(6,*) 'rmax1 < rmin1  iseed = ',iseed
        write(6,*) 'r(k),a1,e1,rmax1,rmin1,m = ',
     &             r(k),a1,e1,rmax1,rmin1,m
        call flush(6)
*
      endif
*
      gmax1 = sqrt(-3.0d0*(rmax1 - rmin1)/gmax1)
      gmin1 = sqrt(3.0d0*(rmax1 - rmin1)/gmin1)
*
      if(rmax1.gt.rmaxz) irm1 = 1
*     
      if(m.gt.n) then
        kmax1 = n
      else
        kmax1 = m + 1
      endif
*
      if((m.eq.1).or.(imi1.eq.3)) kmax1 = maxk
*
C       print*,ipoin1,k,e1,a1,maxk,rmaxz,gmax1,gmin1,istar,rmin1,
C      &     rmax1,imi1,imk1,kmax1,kmin1,irm1
c      stop
*       second star
*
   30 continue
*
      xgmin(im1) = gmin1
      xgmax(im1) = gmax1
*
      if(ipoin2.eq.1) go to 60
*
      m = mmin(k+1,e2,a2,n)
      if (m.gt.0) goto 43
      i = -1
   40 i = i + 1
      m = k + 1 - i
*
      if(m.eq.0) then
        b1 = 0.0d0
        c1 = u(1)
        go to 45
      endif
*
      q1 = 2.0d0*(e2 - u(m)) - a2**2/r(m)**2
      if(q1.ge.0.0d0) go to 40
*
      if(m.eq.k+1) then
        imk2 = 1
        rmin2 = 0.9999d0*r(k+1)
        write(6,*) '< r(k+1) k+1,im2,r,u,e2,vr,vt=',
     &  k+1,im2,r(k+1),u(k+1),e2,vr(im2),vt(im2)
	call flush(6)
        kmin2 = 1
        go to 47
      endif
*
      if(m.eq.n) then
        b1 = u(m)*r(m)
        c1 = 0.0d0
        go to 45
      endif
*
 43   continue
      dz = 1.0d0/r(m+1) - 1.0d0/r(m)
      b1 = (u(m+1) - u(m))/dz
      c1 = (u(m)/r(m+1) - u(m+1)/r(m))/dz
   45 rmin2 = (-b1 + sqrt(b1*b1 - 2.0d0*a2*a2*(c1 - e2)))/a2/a2
      rmin2 = 1.0d0/rmin2
*
      if(rmin2.gt.r(k+1)) then
*
        xrr = 1.d0 - rmin2/r(k+1)
        xrr = abs(xrr)
*
        if(xrr.gt.1.0d-6) then
          write(6,*) '    rmin2  >  r(k+1)  iseed = ',iseed
          write(6,*) 'm,k1,rm2,r(k+1),r(m),r(m+1),u(m),u(m+1),e2,a2 =',
     &                m,k+1,rmin2,r(k+1),r(m),r(m+1),u(m),u(m+1),e2,a2
          call flush(6)
*
          rmin2 = 0.9999d0*r(k+1)
          imi2 = 1
        else
          rmin2 = r(k+1)
          imi2 = 1
        endif
*
      endif
*
      gmin2 = 2.0d0*a2*a2/rmin2**3 + 2.0d0*b1/rmin2**2
*
      if(gmin2.lt.0.0d0) then
*         
        write(6,*) '    gmin2 < 0.0   iseed,gmin2 = ',iseed,gmin2
        write(6,*) 'm,k1,rmin2,r(k+1),r(m),r(m+1),u(m),u(m+1),e2,a2 =', 
     &              m,k+1,rmin2,r(k+1),r(m),r(m+1),u(m),u(m+1),e2,a2
        call flush(6)
*         
        gmin2 = 1.d10   
      endif
*
      if(m.le.1) then
        kmin2 = 1
      else
        kmin2 = m - 1
      endif
*
      if((m.eq.n).or.(imi2.eq.1)) kmin2 = 1
*
 47   continue
*
      m = mmax(k+1,e2,a2,n)
      if (m.gt.0) go to 53

      i = -1
*
   50 i = i + 1
      m = k + 1 + i
*
      if(m.gt.n) then
        b1 = u(m-1)*r(m-1)
        c1 = 0.0d0
        go to 55
      endif
*
      q1 = 2.0d0*(e2 - u(m)) - a2**2/r(m)**2
      if(q1.ge.0.0d0) go to 50
*
      if(m.eq.k+1) then
        if(imk2.eq.1) then
          rnew2 = r(k+1)
          imk2 = 2
          write(6,*) '-- O -- k+1,im2,r,u,e2=',k+1,im2,r(k+1),u(k+1),e2
	  call flush(6)
          rmin2 = rnew2
          rmax2 = rnew2
          go to 60
        endif

        imk2 = 1
        rmax2 = 1.0001d0*r(k+1)
        write(6,*) '> r(k+1) k+1,im2,r,u,e2 = ',k+1,im2,r(k+1),u(k+1),e2
	call flush(6)
        kmax2 = maxk
        go to 60
      endif
*
 53   continue
      dz = 1.0d0/r(m) - 1.0d0/r(m-1)
      b1 = (u(m) - u(m-1))/dz
      c1 = (u(m-1)/r(m) - u(m)/r(m-1))/dz
   55 rmax2 = 0.5d0*(-b1 + sqrt(b1*b1 - 2.0d0*a2*a2*(c1 - e2)))
     &        /(c1 - e2)
*
      if(rmax2.lt.r(k+1)) then
*
        if(imi2.eq.1) then
*         
          write(6,*) ' circular orbit   iseed,k+1 =',iseed,k+1
	  call flush(6)
          imi2 = 2
          rnew2 = r(k+1)
          rmin2 = r(k+1)
          rmax2 = r(k+1)
          go to 60
        endif
*      
        xrr = 1.0d0 - rmax2/r(k+1)
        xrr = abs(xrr)
*
        if(xrr.gt.1.0d-6) then     
          write(6,*) 'rmax2 < r(k+1), rmax2=1.0001*r(k+1),iseed= ',iseed
          write(6,*) 'm,k1,rmax2,r(k+1),r(m),r(m+1),u(m),u(m+1),e2,a2=',
     &              m,k+1,rmax2,r(k+1),r(m),r(m+1),u(m),u(m+1),e2,a2            
          call flush(6)
*                   
          rmax2 = 1.0001d0*r(k+1)
          imi2 = 3
        else
          rmax2 = r(k+1)
          imi2 = 3
        endif
*
      endif
*
      gmax2 = 2.0d0*a2*a2/rmax2**3 + 2.0d0*b1/rmax2**2
*
      if(gmax2.gt.0.0d0) then
*                   
        write(6,*) '    gmax2  >  0.0   iseed,gmax2 = ',iseed,gmax2
        write(6,*) 'm,k1,rmax2,r(k+1),r(m),r(m-1),u(m),u(m-1),e2,a2 =',
     &              m,k+1,rmax2,r(k+1),r(m),r(m-1),u(m),u(m-1),e2,a2            
        call flush(6)
*                   
        gmax2 = -1.d10
      endif
*
      if(rmax2.lt.rmin2) then
*
        write(6,*) 'rmax2 < rmin2  iseed = ',iseed
        write(6,*) 'r(k+1),a2,e2,rmax2,rmin2,m = ',
     &             r(k+1),a2,e2,rmax2,rmin2,m
        call flush(6)
      endif
*
      gmax2 = sqrt(-3.0d0*(rmax2 - rmin2)/gmax2)
      gmin2 = sqrt(3.0d0*(rmax2 - rmin2)/gmin2)
*
      if(m.gt.n) then
        kmax2 = n
      else
        kmax2 = m + 1
      endif
*
      if((m.eq.2).or.(imi2.eq.3)) kmax2 = maxk
*
      if(rmax2.gt.rmaxz) irm2 = 1
*
   60 continue
*
      xgmin(im2) = gmin2
      xgmax(im2) = gmax2
*
C       print*,ipoin2,k,e2,a2,maxk,rmaxz,gmax2,gmin2,istar,rmin2,
C      &     rmax2,imi2,imk2,kmax2,kmin2,irm2
c      stop

*       determination of the new position
*
*       first star
*
      if(ipoin1.eq.1) go to 80
      xmin(im1) = rmin1
      xmax(im1) = rmax1
      xmaxz(im1) = rmaxz
      if(imi1.eq.2) go to 80
      if(imk1.eq.2) go to 80
*
      if(kmax1-kmin1.lt.2) then
        kmax1 = n
        kmin1 = 1
      endif  
*
      ipp = 0
      qmax = 1.9d0
      xcut = 0.01d0
      ycut = 1.d0 - xcut
      rles = rmin1
      rtot = rmax1
      if(irm1.eq.1) then
        rtot = rmaxz
        ycut = 1.0d0
      endif
*
   70 s = ran2(irun)
*
      if(s.lt.xcut) s = xcut
      if(s.gt.ycut) s = ycut
*
      s = -1.0d0 + 2.0d0*s
*
      rnew1 = 0.5d0*(rmin1 + rmax1) + 
     &                    0.25d0*(rmax1 - rmin1)*(3.0d0*s - s**3)
*
      ipp = ipp + 1
*
      if(ipp.lt.10000) then
*
        if((rnew1.gt.rtot).or.(rnew1.lt.rmglob)) go to 70
        if((k.eq.1).and.(rnew1.lt.r(1))) go to 70
        if(ikind(im1).lt.0) then
          print*,'newpos-1: i,im1,ikind = ',i,im1,ikind(im1)
        endif
        if(ikind(im1).eq.2.and.rnew1.lt.r(1)) go to 70
c        if(rnew1.lt.r(1)) go to 70
*
        q0 = qmax*ran2(irun)*dmax1(gmin1,gmax1)
c        print*,'q0,s ',q0,s
        ipo = 0
        rnn1 = rnew1
        call potent(ipo,kmin1,kmax1,rnn1,u1)
        v1 = 2.0d0*(e1 - u1) - a1*a1/rnew1/rnew1
*
      else
*
        write(6,*) 'ipp>10000 isee,time,rmz = ',iseed,time,rmaxz
        write(6,*) 'k,im1,rx1,ri1,v1,q0,e1,u1,a1,rnew1 = ',
     &              k,im1,rmax1,rmin1,v1,q0,e1,u1,a1,rnew1
        call flush(6)
        rnew1 = r(k)
        vrp(im1) = vr(im1)
        vtp(im1) = vt(im1)
        go to 80
*       
      endif
*
c
c      open(33,file='choose.dat',access='append')
c      write(33,1234) time,ipp,k,im1,xcut,ycut,s,rmin1,rnew1,
c     &               rmax1,rmaxz,e1,a1,u1,v1
c 1234 format(1x,'-1- ',1pe12.4,3i6,1p13e12.4)
c      close(33)
c
      if(v1.le.0.0d0) go to 70
      v1 = sqrt(v1)
      drds = 0.75d0*(rmax1 - rmin1)*(1.0d0 - s*s)
      q1 = drds/v1
*
      if(q1.gt.qmax*dmax1(gmin1,gmax1)) then
c        write(6,*) 'q1 > q  iseed,q1,gmin1,gmax1=',iseed,q1,gmin1,gmax1
        qmax = 2.0d0*qmax
        go to 70
      endif
*
c
c      open(33,file='choose.dat',access='append')
c      write(33,1235) time,ipp,k,im1,qmax,gmin1,gmax1,v1,drds,q0,q1
c 1235 format(1x,'-1a- ',1pe12.4,3i6,1p7e12.4)
c      close(33)
c      
      if(q0.gt.q1) go to 70
*
*       check for a nearly radial orbit if the new position is not
*       pick up too close to the apocentre
*
      q1 = a1/rnew1
      vrp(im1) = v1
      vtp(im1) = q1
*
*
   80 continue
*
*       second star
*
      if(ipoin2.eq.1) go to 100
      xmin(im2) = rmin2
      xmax(im2) = rmax2
      xmaxz(im2) = rmaxz
      if(imi2.eq.2) go to 100
      if(imk2.eq.2) go to 100
*
      if(kmax2-kmin2.lt.2) then
        kmax2 = n
        kmin2 = 1
      endif
*
      ipp = 0
      qmax = 1.9d0
      xcut = 0.01d0
      ycut = 1.0d0 - xcut
      rles = rmin2
      rtot = rmax2
      if(irm2.eq.1) then
        rtot = rmaxz
        ytot = 1.0d0
      endif
*
   90 s = ran2(irun)
*
      if(s.lt.xcut) s = xcut
      if(s.gt.ycut) s = ycut
*
      s = -1.0d0 + 2.0d0*s
*
      rnew2 = 0.5d0*(rmin2 + rmax2) + 
     &                       0.25d0*(rmax2 - rmin2)*(3.0d0*s - s**3)
*
      ipp = ipp + 1
*
      if(ipp.lt.10000) then
*
        if((rnew2.gt.rtot).or.(rnew2.lt.rmglob)) go to 90
        if((k+1.eq.2).and.(rnew2.lt.r(1))) go to 90
        if(ikind(im2).lt.0) then
          print*,'newpos-2: i,im2,ikind = ',i,im2,ikind(im2)
        endif
        if(ikind(im2).eq.2.and.rnew2.lt.r(1)) go to 90
c        if(rnew2.lt.r(1)) go to 90
*
        q0 = qmax*ran2(irun)*dmax1(gmin2,gmax2)
        ipo = 0
        rnn2 = rnew2
        call potent(ipo,kmin2,kmax2,rnn2,u2)
        v1 = 2.0d0*(e2 - u2) - a2*a2/rnew2/rnew2
*
      else
*
        write(6,*) 'ipp>10000 isee,t,rmz = ',iseed,time,rmaxz
        write(6,*) 'k+1,im2,rx2,ri2,v1,q0,e2,u2,a2,rnew2 = ',
     &              k+1,im2,rmax2,rmin2,v1,q0,e2,u2,a2,rnew2
        call flush(6)
        rnew2 = r(k+1)
        vrp(im2) = vr(im2)
        vtp(im2) = vt(im2)
        go to 100
*       
      endif
*
c
c      open(33,file='choose.dat',access='append')
c      write(33,1236) time,ipp,k+1,im2,xcut,ycut,s,rmin2,rnew2,
c     &               rmax2,rmaxz,e2,a2,u2,v1
c 1236 format(1x,'-2- ',1pe12.4,3i6,1p13e12.4)
c      close(33)
c
      if(v1.le.0.0d0) go to 90
      v1 = sqrt(v1)
      drds = 0.75d0*(rmax2 - rmin2)*(1.0d0 - s*s)
      q1 = drds/v1
*
      if(q1.gt.qmax*dmax1(gmin2,gmax2)) then
c        write(6,*) 'q1 > q  iseed,q1,gmin2,gmax2 =',iseed,q1,gmin2,gmax2
        qmax = 2.0d0*qmax
        go to 90
      endif
*
c
c      open(33,file='choose.dat',access='append')
c      write(33,1237) time,ipp,k+1,im2,qmax,gmin2,gmax2,v1,drds,q0,q1
c 1237 format(1x,'-2a- ',1pe12.4,3i6,1p7e12.4)
c      close(33)
c
      if(q0.gt.q1) go to 90
*
*       check for a nearly radial orbit if the new position is not 
*       pick up too close to the apocentre
*
      q1 = a2/rnew2
      vrp(im2) = v1
      vtp(im2) = q1
*
*
  100 continue
*
*
      rn1 = rnew1
      rn2 = rnew2
c      print*,rnew1,rnew2
c      stop
      return
*
      end
*
*
*
*
c      subroutine newton(x,r,r2,r1)
*
*       newton iteration for finding minimum and maximum values 
*       -------------------------------------------------------
*       of random number s
*       ------------------
*
*
c      real*8 x,r,r2,r1,ri,f,fp,xs,xs1,dx
*
c      integer iter
*
*
c      iter = 0
c      xs = 0.5d0
*
c      ri = (r - r1)/(r2 - r1)
*
c 10   continue
*
c      f = (3.0d0 - 2.0d0*xs)*xs*xs - ri
c      fp = 6.0d0*xs*(1.0d0 - xs)
*
c      xs1 = xs - f/fp
*
c      dx = abs(xs1 - xs)
c      dx = dx/xs
c      xs = xs1
*
c      iter = iter + 1
c      if(iter.gt.1000) stop  'iter > 1000'
*
c      if(dx.gt.1.0d-6) go to 10
*
c      x = xs
*
c      return
c      end
*
*
*
      function mmin(k,e,a,n)
cThis function returns the m-value for rmin in simple case, otherwise zero
      include 'common.h'
      double precision e,a,qk,q1,qstar,q2
      integer k,n,k1,k2,kstar,nit,mmin
      q1 = 2.0d0*(e - u(1)) - a**2/r(1)**2
      qk = 2.0d0*(e - u(k)) - a**2/r(k)**2
      if(q1.lt.0.d0.and.qk.gt.0.d0.and.k.gt.1.and.k.le.n) then
         k1 = 1
         k2 = k
         q2 = qk
         nit = 1
 10      continue
         if (k2.eq.k1 + 1) then
            mmin = k1
            return
         else
            nit = nit + 1
            if (nit.gt.1000) then
               print*,'nit>1000',k,e,a,n,k1,k2,q1,q2
               mmin = 0
               return
            endif
            kstar = (k1+k2)/2
            qstar = 2.0d0*(e - u(kstar)) - a**2/r(kstar)**2
            if (qstar.lt.0.d0) then
               k1 = kstar
               q1 = qstar
               go to 10
            else if (qstar.gt.0.d0) then
               k2 = kstar
               q2 = qstar
               go to 10
            else
               mmin = 0
               print*,'qstar = 0',k,e,a,n,k1,k2,q1,q2
               return
            endif
         endif
      else
         print*,'non-standard',k,e,a,n,q1,qk
         mmin = 0
         return
      endif
      end

      function mmax(k,e,a,n)
cThis function returns the m-value for rmax in simple case, otherwise zero
cNote that it is the index of the star above the maximu
      include 'common.h'
      double precision e,a,qk,qn,qstar,q1,q2
      integer k,n,k1,k2,nit,kstar,mmax
      q1 = 0.d0
      q2 = 0.d0
      qk = 2.0d0*(e - u(k)) - a**2/r(k)**2
      qn = 2.0d0*(e - u(n)) - a**2/r(n)**2
      if(qn.lt.0.d0.and.qk.gt.0.d0.and.k.ge.1.and.k.lt.n) then
         k1 = k
         k2 = n
         q1 = qk
         q2 = qn
         nit = 1
 10      continue
         if (k2.eq.k1 + 1) then
            mmax = k2
            return
         else
            nit = nit + 1
            if (nit.gt.1000) then
               print*,'nit>1000',k,e,a,n,k1,k2,q1,q2
               mmax = 0
               return
            endif
            kstar = (k1+k2)/2
            qstar = 2.0d0*(e - u(kstar)) - a**2/r(kstar)**2
            if (qstar.gt.0.d0) then
               k1 = kstar
               q1 = qstar
               go to 10
            else if (qstar.lt.0.d0) then
               k2 = kstar
               q2 = qstar
               go to 10
            else
               mmax = 0
               print*,'qstar = 0',k,e,a,n,k1,k2,q1,q2
               return
            endif
         endif
      else
         print*,'non-standard',k,e,a,n,q1,qk
         mmax = 0
         return
      endif
      end

