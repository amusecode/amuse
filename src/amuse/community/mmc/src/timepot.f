*
*
      subroutine timpot(lmax)
*
*
*       determine changes to the system due to changes of mass distribution
*       -------------------------------------------------------------------
*
*     timpot is called before interaction sequence (intb3b3 - intb3f - formb3)
*     so ikind cannot be negative
*
*
      include 'common.h'
*
c      real*8 uoo,uon,uno,unn,vold,vnew,dd,v2t,vrp,vtp
      real*8 uoo,uon,uno,unn,vold,vnew,dd,v2t
*  
      integer lmax,i,k1,k2,im,imo,ipo,iex,im1
*
c      common /timep/ vrp(nmax),vtp(nmax)
*
      do 10 i=1,lmax
         im = iname(i)
         imo = inameo(i)
*
*     determine old potentials for old and new positions
*
         uoo = uo(imo)
         ipo = 1
         k1 = 1
         k2 = lmax
         call potent(ipo,k1,k2,r(i),uon)
*
*     determine new potentials for old and new positions
*
         ipo = 0
         k1 = 1
         k2 = lmax
         call potent(ipo,k1,k2,ro(imo),uno)
         unn = u(i)
*
*     determine old and new velocities
*
         vold = vro(im)*vro(im) + vto(im)*vto(im)
         vnew = vold + uoo - uon - unn + uno
* 
*     determine new radial and tangential velocities
*
         vt(im) = vto(im)*ro(imo)/r(i)
         v2t = vt(im)*vt(im)
         dd = vnew - v2t
         if(body(im).lt.0.d0) then
           write(6,*) ' m < 0  i,im,imo,m = ',i,im,imo,body(im)
           stop ' mass < 0 '
         endif
*
         if(dd.gt.0.0d0) then
           vr(im) = sqrt(dd)
           vrr(im) = 0.0d0
         else   
*
           if(vnew.lt.0.0d0) then
             vr(im) = 0.0d0
             vrr(im) = 0.5d0*dd*body(im)
           else
             dd = xmax(im)/r(i) - 1.0d0
             v2t = r(i)/xmin(im) - 1.0d0
             if((dd.lt.0.03d0).or.(v2t.lt.0.03d0)) then
               vr(im) = 0.0d0
               vt(im) = sqrt(vnew)
               vrr(im) = 0.0d0
             else
               dd = vrp(im)**2 + vtp(im)**2
               dd = sqrt(dd/vnew)
               vr(im) = 0.0d0
               vt(im) = sqrt(vnew)
               vr(im) = vrp(im)/dd
               vt(im) = vtp(im)/dd
               vrr(im)= 0.0d0
             endif
           endif
         endif
*
*     for the nminzo stars the velocities are artificialy isotropised
*
c 20      continue
c  20     iex = 3*nminzo
c  20     iex = 0.1*nt
         iex = -10
         if(i.le.iex) then
         im1 = iname(i)
         if(ikind(im1).lt.0) then
           print*,'timpot: i,im1,ikind = ',i,im1,ikind(im1)
         endif
           if(ikind(im1).eq.2) go to 10
           dd = (vr(im)**2 + vt(im)**2)/3.0d0
           vr(im) = sqrt(dd)
           vt(im) = sqrt(2.0d0*dd)
         endif
*
 10   continue
*
*
      return
*
      end
*
*
*
