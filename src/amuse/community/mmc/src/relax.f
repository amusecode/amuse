      subroutine relax(nup,lmin,lmax,dt)
*
*
*
*       determine relaxation process for stars between lmin and lamx
*       ------------------------------------------------------------
*       with time-step dt
*       -----------------
*
*
      include 'common.h'
*
      real*8 sm1,sm2,w,den,a12,sin2b2,dt1,beta,dt,dt0,rx1,rx2,ek,ek0,
     &       ep,ep1,ep2,voo,dt2
*
      real*4 ran2
*
      integer i,im1,im2,lmin,lmax,nup,nc2,ibound,m,k,ibx
*
*
      ibound = 0
      ep = 0.0d0
      ek = 0.0d0
      ep1 = 0.0d0
      ep2 = 0.0d0
      dt0 = dt
      nc2 = nc/2
      if(nc2.gt.nminzo) nc2 = nminzo
*
*      determine deflection angle and new velocities
*
      do 10 i=lmin,lmax,2
*
         im1 = iname(i)
*
         vro(im1) = vr(im1)
         vto(im1) = vt(im1)
         ro(i) = r(i)
         rn(i) = r(i)
         if(i+1.gt.lmax) go to 10
*
         im2 = iname(i+1)
         sm1 = body(im1)
         sm2 = body(im2)
cAdded DCH 2/8/6
         if (sm1*sm2.eq.0.d0) go to 10
*
*       do not relax pair of stars which apocentre distances are 
*       smaller then actual position of the next star
*
c        if(i.lt.nc2) then
        if(i.lt.0) then
          if(xmax(im1).eq.0.0d0) go to 15
          rx1 = xmax(im1)
          rx2 = xmax(im2)
          if(rx1.gt.xmaxz(im1)) rx1 = xmaxz(im1)
          if(rx2.gt.xmaxz(im2)) rx2 = xmaxz(im2)
          if((rx1.lt.r(i+2)).and.(rx2.lt.r(i+2))) go to 30
        endif
*
 15     continue
*
*       determine relative velocity of stars i and i+1
*
         if(ran2(irun).lt.0.5) vr(im1) = -vr(im1)
         if(ran2(irun).lt.0.5) vr(im2) = -vr(im2)
*
         call relvel(i,0,w)
*
*       determine density in the vicinity of stars i and i+1
*
         call densi(i,den) 
*
*       determine the deflection angle
*
         a12 = den*(sm1 + sm2)**2/w**3*float(nt)
         sin2b2 = a12*dt0
* 
*       if beta greater than PI form a sequence of interaction
*
         dt1 = dt0
   20    continue
         if(sin2b2.gt.1) then
*
           dt2 = 0.5d0/a12
c           dt2 = 1.0d0/a12
           dt1 = dt1 - dt2
*
*       find new velocities after interaction with deflection angle = PI/2
*
           beta = 0.5d0*pi
c           beta = pi
           call newvel(i,beta,den,dt2)
*
           vr(im1) = vrn1
           vr(im2) = vrn2
           vt(im1) = vtn1
           vt(im2) = vtn2
*
           call relvel(i,0,w)
*
*       determine new deflection angle
*
           a12 = den*(sm1 + sm2)**2/w**3*float(nt)
           sin2b2 = a12*dt1
c           dt1 = dt2
           go to 20
         endif
*
         beta = 2.0d0*atan(sqrt(sin2b2/(1.0d0-sin2b2)))
*
*       determine velocities after interaction
*
         call newvel(i,beta,den,dt1)
*
         vr(im1) = vrn1
         vr(im2) = vrn2
         vt(im1) = vtn1
         vt(im2) = vtn2
*
*       compute the self-energy of up to nminzo innermost stars
*       if it is negative devide total kinetic energy of these stars
*       on two equal parts and isotropise the velocities
*
 30     continue
        ibound = 1
        if(ibound.eq.0) then
          if(i.le.nminzo) then
            ek0 = 0.5d0*(sm1*(vr(im1)**2 + vt(im1)**2) +
     &                   sm2*(vr(im2)**2 + vt(im2)**2))
            ek = ek + ek0
            ep1 = ep1 + 0.5d0*(sm1**2/r(i) + sm2**2/r(i+1))
            do 40 k=2,i+1
            do 50 m=k,i+1
               ep2 = ep2 + body(iname(k-1))*body(iname(m))/r(m)
 50         continue
 40         continue
            ep = ep1 + ep2
            ep2 = 0.0d0
            if(ep.lt.ek) then
              ibound = 1
              if((i+1).eq.2) go to 10
              ek0 = (ek - ek0)/float(i-1)
              do 60 k = 1,i-1
                 im1 = iname(k)
cChanged by DCH 2/8/6
                 if (body(im1).ne.0.d0) then
                    voo = 2.0d0*ek0/3.0d0/body(im1)
                    vt(im1) = sqrt(2.0d0*voo)
                    vr(im1) = sqrt(voo)
                 endif
 60           continue
*
              write(6,*) 'bound num,ep,ek,lmi,lma =',i-1,ep,ek,lmin,lmax
*
              go to 10
            endif
          endif
        endif
*
 10   continue  
*
*
      do 70 i=lmin,lmax,2
*
*       define velocities of stars after interaction as old velocities
*
         if(i+1.gt.lmax) go to 70
*
         im1 = iname(i)
         im2 = iname(i+1)         
*
         vro(im1) = vr(im1) 
         vro(im2) = vr(im2)
         vto(im1) = vt(im1)
         vto(im2) = vt(im2)
*
*       determine new positions of two interacting stars
*
c
c         rrroo1 = r(i)
c         rrroo2 = r(i+1)

         call newpos(i,nup)
*
*       define positions of stars before interaction as old position
*
         ro(i) = r(i)
         ro(i+1) = r(i+1)
         rn(i) = rn1
         rn(i+1) = rn2
*
*        set new position for binaries (binary array)
*
         if(ikind(im1).ge.2) then
c           ibx = im1 - nss0
c           if(ibx.lt.0) then
           ibx = nwhich(im1)
c           endif
           bin(ibx,8) = rn1
           bin(ibx,6) = time + dt
         endif
*
         if(ikind(im2).ge.2) then
c           ibx = im2 - nss0
c           if(ibx.lt.0) then
           ibx = nwhich(im2)
c           endif
           bin(ibx,8) = rn2
           bin(ibx,6) = time + dt
         endif
*                                                                                    *
c
c         if(rn1.gt.rtid.or.rn2.gt.rtid) then
c      print*,'rrr: im1,im2,i,i+1,ro1,ro2,rn1,rn2,2(rmax,rmin,rmaxz)=',
c     &im1,im2,i,i+1,rrroo1,rrroo2,rn1,rn2,xmin(im1),xmax(im1),
c     &xmaxz(im1),xmin(im2),xmax(im2),xmaxz(im2)
c         endif
*
 70   continue
*
*
      return
*
      end
*
*
*
*
