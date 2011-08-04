      function polardens(r,ctheta,l)
      parameter(pi=3.1415926535)
c returns density at r,theta, multiplied by the lth harmonic
c (must integrate this over theta to get the harmonic coeff)
      common/legendre/ plcon(0:40)
      real psi

      z=r*ctheta
      s=r*sqrt(1.0 - ctheta*ctheta)
      polardens=dens(s,z)*plcon(l)*plgndr1(l,ctheta)

      return
      end
