      subroutine effectivescalelength(aeff)

      include 'commonblocks'
      
      
      sig2 = -psi0
      vb2 = v0bulge*v0bulge
      v02 = v0*v0
      if(v0bulge.ne.0.) then
         coef1 = 4./v0bulge/abulge/abulge*((vb2/sig2)**2.5)
      else
         coef1 = 0.
      endif
      coef2 = 1./v0/a/a*(((v02+2.*vb2)/sig2)**2.5)	  
      coef3 = sig2**1.5
      aeff = sqrt(4.*sig2/coef3/(coef1+coef2));
      return
      end

      subroutine dfcorrectionhalo(e,fac)

      include 'commonblocks'
      
      sig2 = -psi0
      call effectivescalelength(aeff)
c      fac = 1 + 3.*pi/4.*
c     +     ((bhmass/aeff/sig2)**1.5)/((1-e/sig2)**3.)

      del = 0.1/aeff
      den = sqrt((1-e/sig2)**2. + del**2.)
      fac = 1 + 3.*pi/4.*
     +     ((bhmass/aeff/sig2)**1.5)/(den**3.)
      return
      end

      subroutine dfcorrectionbulge(e,fac)

      include 'commonblocks'

      sig2 = -psi0
      call effectivescalelength(aeff)
      fac = 1 + 3.*pi/4.*
     +     ((bhmass/aeff/sig2)**1.5)/((1-e/sig2)**3.)
      return
      end

      function potbh(R,z)

      include 'commonblocks'

      potbh = - bhmass/sqrt(R*R + z*z)

      return
      end

      
