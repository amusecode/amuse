
      function polarbulgedens(r,ctheta,l)

      common/legendre/ plcon(0:40)
c returns density at r,theta, multiplied by the lth harmonic
c (must integrate this over theta to get the harmonic coeff)
      z=r*ctheta
      s=r*sqrt(1.0 - ctheta*ctheta)
      polarbulgedens=bulgedenspsi(pot(s,z))*
     &               plcon(l)*plgndr1(l,ctheta)
      return
      end
