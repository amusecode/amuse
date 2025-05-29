
      function appdiskpot(s,z)
      parameter (pi=3.1415926535)
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea

c  The potential is obtained by integrating the density of a sech**2
c  disk up vertically, ignoring radial potential gradients.
c  In the resulting expression the cylindrical r is then replaced with
c  spherical radius, to make the potential zero at infinity 
c  (otherwise the harmonics series do not have a simple asymptotic 
c  form at infinity).
c      appdiskpot = 0
c      return
      r=sqrt(s*s+z*z)
c  erfc radial truncation factor
      t=sqrt(0.5d0)*(r-outdisk)/drtrunc
      if (t.lt.-4.0) then
        eerfc=1.
      elseif(t.gt.4.0) then
        appdiskpot=0.
        return
      else
        eerfc=0.5*erfc(t)
        endif

      f=diskconst*exp(-r/rdisk)
      zz=abs(z/zdisk)
      appdiskpot=4*pi*f*zdisk**2*(zz+log(0.5*(1+exp(-2*zz))))*eerfc
      return
      end       
