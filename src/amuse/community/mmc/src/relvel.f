      subroutine relvel(kn,kk,w)
*
*
*       determine relative velocity of two interacting stars
*       ----------------------------------------------------
*
*
      include 'common.h'
*
*
      real*8 w,vr1,vr2,vtx1,vtx2,vty1,vty2,wx,wy,wz,vx,vy,vz,
     &       sm1, sm2, sm12
*
      real*4 ran2
*
      integer k,kk,kn,im1,im2
*
*
      k = abs(kn)
      im1 = iname(k)
*
      if(kk.eq.0) then
        im2 = iname(k+1)
      else
        im2 = iname(kk)
      endif
*
      vr1 = vr(im1)
      vr2 = vr(im2)
      sm1 = body(im1)
      sm2 = body(im2)
      sm12 = sm1 + sm2
*
      fi = twopi*ran2(irun)
      sinfi = sin(fi)
      cosfi = cos(fi)
*
      vtx1 = vt(im1)
      vtx2 = vt(im2)*cosfi
      vty1 = 0.0d0
      vty2 = vt(im2)*sinfi
*
      wx = vtx2 - vtx1
      wy = vty2 - vty1
      wz = vr2 - vr1
      w = wx*wx + wy*wy + wz*wz
      w = sqrt(w)
*
*     compute ceter of mass velocity components for mergrer
*     k < 0 means merger 
*
      if(kn.lt.0) then
        vx = (sm1*vtx1 + sm2*vtx2)/sm12
        vy = (sm1*vty1 + sm2*vty2)/sm12
        vz = (sm1*vr1 + sm2*vr2)/sm12
        vr(im1) = vz
        vt(im1) = sqrt(vx*vx + vy*vy)
      endif
*
      return
*
      end
*
*
*
*
      
