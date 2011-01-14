      subroutine zone
*
*
*       compute devision of the system on spherical zones and grup
*       ----------------------------------------------------------
*       them in super-zones. Compute individual time steps for
*       ------------------------------------------------------
*       each super-zone.
*       ----------------
*
*
*
      include 'common.h'      
*
*
      real*8 smz,w2,dt,beta,xld,az,tref,tccc,xxl,bmin1,bmax1,vv,rr
*
      integer n,i,nzon,num,nco,ibound,nud,lo,ln,ld,k,im,itwo
*
      dimension az(nzonma),vv(nzonma),rr(nzonma)
*
*
      n = nt
*
*       determine number of zones and stars in each zone
*
      if (nt/nz0.gt.nzonma) then
         print*,'zone: nt/nz0 > nzonma'
         stop
      endif
      if(time.eq.0.0d0) then
        nzon0 = nt/nz0
        izon = 0
        nzon = nzon0
      else
        nzon = n/nz0
        if(izon.eq.1) go to 10
*
        if(nzon.lt.0.9*nzon0) then
          nz0 = nz0 - 2
          izon = 0
          nzon = n/nz0
*
          if(nz0.lt.nminzo) then
            nz0 = nminzo
            nzon = n/nz0
            izon = 1
          endif
*
        endif
*
      endif
*
   10 continue
*
*       determine the name of the last star in each zone
*
      i = 0
   20 i = i + 1
      if(i.gt.nzon-1) go to 30
      nzs(i) = nz0*i
      go to 20
*
   30 nzs(i) = n 
*
*
*       check if in the core is at least nzonc zones
*
      i = 0
   40 i =i + 1
      if(nzs(i).lt.nc) go to 40
*
      if(i.ge.nzonc) go to 90
*
      ibound = i
      nco = nzs(ibound)
      num = nco/nzonc
*
        if(num.gt.nminzo) then
*
          nud = nzonc - ibound
          k = mod(num,2)
*
            if(k.gt.0) then
              num = num - 1
            endif
*
          do 50 i=nzon,ibound+1,-1
   50        nzs(i+nud) = nzs(i)
*
          do 60 i=1,nzonc-1
   60        nzs(i) = num*i
*
          nzs(nzonc) = nco
          nzon = nzon + nud
*
        else
* 
          num = nzs(ibound)/nminzo
          nud = num - ibound
*
          do 70 i=nzon,ibound+1,-1
   70        nzs(i+nud)=nzs(i)
*
          do 80 i=1,num-1
   80        nzs(i)= nminzo*i
*
          nzs(num) = nco
          nzon = nzon + nud
*
        endif
*     
*
   90 continue
*
*       determine for each zone coefficient in formulae for sin^2(beta/2)
*       Stodolkiewicz 1982, A.A. Vol.32, 69, eq.(20)
*
      lo = 1   
*
      do 200 i=1,nzon
         ln = nzs(i)
         ld = ln + 1 - lo
         smz = 0.0d0
         w2 = 0.0d0
*
      do 210 k=lo,ln
         im=iname(k)
         smz = smz + body(im)
  210    w2 = w2 + body(im)*(vr(im)*vr(im) + vt(im)*vt(im))
*
         xld = float(ld)
         smz = smz/xld
         w2 = 2.0d0*w2/smz/xld
         vv(i) = sqrt(0.5d0*w2)
         rr(i) = 0.5d0*(r(lo) + r(ln))
*
         az(i) = 6.0d0*nt*xld*smz*smz/w2**1.5/(r(ln)**3 - r(lo)**3)
*
         if(i.eq.2) then
           if(az(1).gt.2.d0*az(2)) az(1) = 1.9d0*az(2)
         endif
*
         if(i.gt.3) then
           if((az(i-2).lt.az(i-1)).and.(az(i).lt.az(i-1))) 
     &                     az(i-1) = 0.5d0*(az(i-2) + az(i))
         endif
*
         if(i.gt.3) then
           if((az(i-2).gt.az(i-1)).and.(az(i).gt.az(i-1))) 
     &                     az(i-1) = 0.5d0*(az(i-2) + az(i))
         endif
*
         lo = ln + 1
* 
  200 continue
*
*       determine zones for which time-step is the same and collect them
*       in a one super-zone.
*       determine time step for each super-zone
*
c      tau = tau0
      xxl = log(gamma*nt)/nt
c      tccc = rc*xxl/vc
      bmin1 = bmin
      bmax1 = bmax
*
  300 k = 0
*
      do 310 i=1,nzon
         itwo = 0
*
  320    itwo = itwo + 1
         if(itwo.gt.ntwo) go to 370
         dt = az(i)*2.0d0**(-itwo)
         beta = dt*tau
*
         if((itwo.eq.1).and.(beta.lt.bmin1)) then
           itwo = 0
           go to 325
         endif
*
         if(beta.gt.bmax1) go to 320
*
  325    continue
*
c         if(itwo.gt.0) then
c           tref = 2.0d0**(-itwo)*tau
c           tccc = rr(i)*xxl/vv(i)
c           if(tref.lt.1.5d0*tccc) then
c           if((k.le.3).and.(tref.lt.1.5d0*tccc)) then
c	      bmin1 = bmax1
c	      bmax1 = 2.0d0*bmax1
c	      open(23,file='crostime.dat',access='append')
c              write(23,*) time,tref,tccc,bmax,tau,bmax1,itwo,k
c	      close(23)
c	      tau = tau0
c	      go to 300
c            endif
c         endif
*
         ln = nzs(i)
         k = k + 1
         nzst(k) = ln
         ltwo(k) = itwo
*
         if(k.eq.1) go to 310
*
         ld = ltwo(k) - ltwo(k-1)
*
         if(ld) 330,340,350
*
  330    if(-ld.eq.1) go to 310
         ltwo(k) = ltwo(k-1) - 1
         go to 310
*
  340    nzst(k-1) = nzst(k)
         k = k - 1
         go to 310
*
  350    continue
* 
         nzst(k-1) = nzst(k)
         k = k - 1
         go to 310
*
  370    tau = tau/2.0d0
         go to 300
*
  310 continue
*
      nsuzon = k
c      nsuzon = ltwo(1)
*
*       print control data
*
c      if(iprint.eq.0) then
        write(6,*) ' nzon = ',nzon,'  nsuzon = ',nsuzon,'  tau = ',tau
*
c        do 400 k=1,nzon
c           bbb = az(k)*tau*2.d0**(-nsuzon+1)
c  400      write(6,*) k,nzs(k),az(k),bbb
*
        do 410 k=1,nsuzon
  410      write(6,*) k,nzst(k),ltwo(k)
*
c      endif
*
*
      return
*
      end
*
*
*
