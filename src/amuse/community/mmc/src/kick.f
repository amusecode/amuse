*
*
      subroutine kick(k,vb)
*
*      compute components of the 'kick' velocity (gained during
*      --------------------------------------------------------
*      interaction between objects) of object and new radial
*      ------------------------------------------------------- 
*      and tangencial velocities
*      -------------------------

      include 'common.h'
*
      real*8 q1,q2,cteta,steta,cosffi,sinffi,vb,b1,c1,dv,zz,xvr,xvt
*
      real*4 ran2
*
      integer im1,k
*
*      used formulae from Stodolkiewicz Acta Astronomica 1986, 36, p19
*
      q1 =ran2(irun) - 0.5d0
      cteta = -2.0d0*q1
      steta = sqrt(1.0d0 - cteta**2)
      im1 = iname(k)
      xvr = vr(im1)
      xvt = vt(im1)
*
      q2 = 2.0d0*pi*ran2(irun)
      sinffi = sin(q2)
      cosffi = cos(q2)
      b1 = 2.0d0*(vr(im1)*cteta + vt(im1)*steta*cosffi)
      c1 = -2.0d0*vb/body(im1)
      dv = 0.5d0*(-b1 + sqrt(b1*b1 - 4.0d0*c1))
      vr(im1) = vr(im1) + dv*cteta
      b1 = vt(im1) + dv*steta*cosffi
      c1 = dv*steta*sinffi
      vt(im1) = sqrt(b1*b1 + c1*c1)
*
      zz = vt(im1)**2 + vr(im1)**2 - xvr**2 - xvt**2
      zz = 0.5d0*body(im1)*zz - vb
*
c      if(vb.gt.1.0d-4) then
c        write(6,*) 'k,im1,r,xvr,xvt,vr,vt,zz =',k,im1,r(k),xvr,xvt,
c     &              vr(im1),vt(im1),zz
c      endif
*
      return
      end
*
*
 
       
