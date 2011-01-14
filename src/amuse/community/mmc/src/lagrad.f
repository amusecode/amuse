      subroutine lagrad
* 
*
*       lagrangian radii & mean radial and tangential velocities in the 
*       ---------------------------------------------------------------
*       lagrangian shells
*       -----------------
*
*
      include 'common.h'
*
      real*8 zm,zm1,zm2,vvr,vvt,zmi,vr2,vt2,r21,r22,rla
*
      integer i,j,n,im
*
*
*       determine the lagrangian radii and mean radial and tangential
*       velocities in the lagrangian shells
*
      n = nt
*
      im = 0
      zm = 0.0d0
      zm2 = 0.0d0
*
      do 10 j=1,nlagra
         vvr = 0.0d0
         vvt = 0.0d0
         zmi = flagr(j)*smt
   20    im = im + 1
         i = iname(im)
         zm = zm + body(i)
         vr2 = vr(i) * vr(i)*body(i)
         vt2 = vt(i) * vt(i)*body(i)
         vvr = vvr + vr2
         vvt = vvt + vt2
*
         if(im.eq.n) then
           vr2 = 0.0d0
           vt2 = 0.0d0
           go to 30
         endif
*
         if (zm.le.zmi) go to 20
   30    zm1 = zm - body(i)
         vvr = vvr - vr2
         vvt = vvt - vt2
         if(im.eq.1) then
           r21 = 0.d0
         else
           r21 = r(im - 1)**3
         endif
         v2rl(j) = vvr/(zm1 - zm2)
         v2tl(j) = vvt/(zm1 - zm2)
         ani(j) = 2.0d0 - v2tl(j)/v2rl(j)
         r22 = r(im)**3
         rla = (zmi - zm) * (r22 - r21)/(zm - zm1) + r22
         rlag(j) = rla**one3
         im = im - 1
         zm = zm1
         zm2  =  zm1
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
*
