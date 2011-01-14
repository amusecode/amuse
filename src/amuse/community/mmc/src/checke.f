      subroutine checke(k,dt,den)
*
*
*       determine if escaping stars can realy escape from the system.
*       -----------------------------------------------------------------
*       Method proposed by D.C. Heggie - 1996.  
*       -----------------------------------------------------------------
*
*
      include 'common.h'
*
*
      real*8 vr1,vr2,vtx1,vtx2,vty1,vty2,beta,wx,wy,
     &       wz,w,wx1,wy1,wz1,wx2,wy2,wz2,sinpsi,cospsi,
     &       sbeta,cbeta,wnx,wny,wnz,vnx1,vnx2,vny1,vny2,
     &       sm12,sm112,sm212,dt,den,lambda,cden,yran,pp,po,pppo,
     &       e1,e2,vrn1n,vrn2n,vtn1n,vtn2n,sin2b2
*
      real*4 ran2
*
      integer k
*
      common /chesc/ sm12,sm112,sm212,w,wx,wy,wz,wx1,wy1,wz1,wx2,
     &               wy2,wz2,vtx1,vtx2,vty1,vty2,vr1,vr2,
     &               sinpsi,cospsi
*
*
*     check if stars escape after relaxation
*
      e1 = u(k) + 0.5d0*(vrn1**2 + vtn1**2)
      e2 = u(k+1) + 0.5d0*(vrn2**2 + vtn2**2)
*
      if((e1.gt.0.0d0).or.(e2.gt.0.0d0)) then
*
*     determine deflection angle beta
*
      cden = 0.5d0/pi
      lambda = cden*pi*den*w*dt*nt/log(gamma*nt)
      yran = ran2(irun) - 1.d0
      yran = abs(yran)
      pp = sqrt(-log(yran)/lambda)
      po = sm12/w/w
      pppo = (pp/po)**2
      sin2b2 = 1.0d0/(1.0d0 + pppo)
      beta = 2.0d0*atan(sqrt(sin2b2/(1.0d0 - sin2b2)))
*
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
      vrn1n = vr1 - sm212*(wnz - wz)
      vrn2n = vr2 + sm112*(wnz - wz)
      vtn1n = sqrt(vnx1*vnx1 + vny1*vny1)
      vtn2n = sqrt(vnx2*vnx2 + vny2*vny2)
*
*     check if any star has positive energy - escapers
*
c      if((e1.gt.0.0d0).or.(e2.gt.0.0d0)) then
        vrn1 = vrn1n   
        vtn1 = vtn1n
        vrn2 = vrn2n
        vtn2 = vtn2n
      endif
*
c      e1 = u(k) + 0.5d0*(vrn1n**2 + vtn1n**2)
c      e2 = u(k+1) + 0.5d0*(vrn2n**2 + vtn2n**2)
*
c      if((e1.gt.0.0d0).or.(e2.gt.0.0d0)) then
c        vrn1 = vrn1n
c        vtn1 = vtn1n
c        vrn2 = vrn2n
c        vtn2 = vtn2n
c      endif
*
      return
*
      end
*
*
*
*
      
