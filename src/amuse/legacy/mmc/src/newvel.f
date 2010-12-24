      subroutine newvel(k,beta,den,dt)
*
*
*       determine new velocities of two interacting objects after encounter
*       -------------------------------------------------------------------
*       M. HENON 1971, Astrophysics and Space Science Vol. 14, 151
*       ----------------------------------------------------------
*
*
      include 'common.h'
*
*
      real*8 vr1,vr2,vtx1,vtx2,vty1,vty2,beta,wx,wy,
     &       wz,w,wx1,wy1,wz1,wx2,wy2,wz2,wp,psi,sinpsi,cospsi,
     &       sbeta,cbeta,wnx,wny,wnz,vnx1,vnx2,vny1,vny2,sm1,sm2,
     &       sm12,sm112,sm212,den,dt
*
      real*4 ran2
*
      integer k,im1,im2
*
      common /chesc/ sm12,sm112,sm212,w,wx,wy,wz,wx1,wy1,wz1,wx2,
     &               wy2,wz2,vtx1,vtx2,vty1,vty2,vr1,vr2,
     &               sinpsi,cospsi
*
*
      im1 = iname(k)
      im2 = iname(k+1)
*
      sm1 = body(im1)
      sm2 = body(im2)
      sm12 = sm1 + sm2
      sm112 = sm1/sm12
      sm212 = sm2/sm12
      vr1 = vr(im1)
      vr2 = vr(im2)
*
      vtx1 = vt(im1)
      vtx2 = vt(im2)*cosfi
      vty1 = 0.0d0
      vty2 = vt(im2)*sinfi
*
      wx = vtx2 - vtx1
      wy = vty2 - vty1
      wz = vr2 - vr1
      wp = wx*wx + wy*wy
      w = wp + wz*wz
      w = sqrt(w)
      wp = sqrt(wp)
*
      wx1 = wy*w/wp
      wy1 = -wx*w/wp
      wz1 = 0.0d0
      wx2 = -wx*wz/wp
      wy2 = -wy*wz/wp
      wz2 = wp
*
      psi = twopi*ran2(irun)
      sinpsi = sin(psi)
      cospsi = cos(psi)
      sbeta = sin(beta)
      cbeta = cos(beta)
*
      wnx = wx*cbeta + wx1*sbeta*cospsi + wx2*sbeta*sinpsi
      wny = wy*cbeta + wy1*sbeta*cospsi + wy2*sbeta*sinpsi
      wnz = wz*cbeta + wz1*sbeta*cospsi + wz2*sbeta*sinpsi
*
      vnx1 = vtx1 - sm212*(wnx - wx)
      vnx2 = vtx2 + sm112*(wnx - wx)
      vny1 = vty1 - sm212*(wny - wy)
      vny2 = vty2 + sm112*(wny - wy)
      vrn1 = vr1 - sm212*(wnz - wz)
      vrn2 = vr2 + sm112*(wnz - wz)
      vtn1 = sqrt(vnx1*vnx1 + vny1*vny1)
      vtn2 = sqrt(vnx2*vnx2 + vny2*vny2)
c
c      write(6,*) 'k,im1,im2,vrn1,vtn1,vrn2,vtn2 = ',k,im1,im2,vrn1,
c     &            vtn1,vrn2,vtn2
   
*
*     check for the true escapers. Method suggested by D.C. Heggie - 1996
*     do not check for tidaly limited clouster
*
cChanged for M67:
c      if(imodel.ne.3) call checke(k,dt,den)
c      if(imodel.ne.3.and.imodel.ne.4) call checke(k,dt,den)
cNext change recommended by mig 25/8/9
      if(imodel.lt.3) call checke(k,dt,den)
*
      return
*
      end
*
*
*
*
      
