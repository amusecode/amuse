      function dpolardiskdens(r,ctheta)
      z=r*ctheta
      s=r*sqrt(1.0 - ctheta*ctheta)
      dpolardiskdens=diskdensestimate(s,z)
      return
      end

      function diskdensestimate(s,z)

      include 'commonblocks'
c
      r=sqrt(s*s+z*z)
      t=sqrt(0.5d0)*(r-outdisk)/drtrunc
      t2=t*t
      if (t.lt.-5.0) then
         eerfc=1.
      elseif (t.lt.5.0) then
         eerfc=0.5*erfc(t)
      else
         eerfc=0
      endif
      if (r.gt.0.) then
         f=diskconst*eerfc*exp(-r/rdisk)
      else
         f=diskconst*eerfc
      endif
      zz=abs(z/zdisk)
      ezz=exp(-zz)
      e2zz=ezz*ezz
      tsech2=(2*ezz/(1+e2zz))**2
      diskdensestimate=f*tsech2

      return
      end

      function halodensity(r)

      include 'commonblocks'
      real nfwdens

      halodensity = nfwdens(r)*eerfc(r)

      return
      end

      function halodensprime(r)

      include 'commonblocks'
      real nfwdens, nfwdensprime

      halodensprime = nfwdens(r)*eerfcprime(r)+nfwdensprime(r)*eerfc(r)
      
      return
      end

      function halodens2prime(r)

      include 'commonblocks'
      real nfwdens, nfwdensprime, nfwdens2prime

      t1 = nfwdens2prime(r)*eerfc(r)
      t2 = 2.*nfwdensprime(r)*eerfcprime(r)
      t3 = nfwdens(r)*eerfc2prime(r)
      halodens2prime = t1 + t2 + t3
      return
      end

      function eerfc(r)

      include 'commonblocks'

      t=sqrt(0.5)*(r-chalo)/drtrunchalo
      t2=t*t
      if (t.lt.-4.) then
         eerfc=1.
      elseif (t.lt.4.) then
         eerfc=0.5*erfc(t)
      else
         eerfc=0
      endif
      
      return
      end

      function eerfcprime(r)

      include 'commonblocks'

      t=sqrt(0.5)*(r-chalo)/drtrunchalo
      t2=t*t

      if(t2.gt.16.) then
         eerfcprime = 0.
      else
         eerfcprime = -0.5*sqrt(2./pi)/drtrunchalo*exp(-t2)
      endif
      return
      end

      function eerfc2prime(r)

      include 'commonblocks'

      t=sqrt(0.5)*(r-chalo)/drtrunchalo
      t2=t*t

      if(t2.gt.16.) then
         eerfc2prime = 0.
      else
         eerfc2prime = 1./sqrt(pi)/drtrunchalo/drtrunchalo*t*exp(-t2)
      endif
C      write(77,*) r,eerfc2prime
      return
      end
