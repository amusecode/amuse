      function appdiskdens(s,z)

      include 'commonblocks'

c  This is the density corresponding to the first-guess disk potential
c  f(r)*erfc((r-outdisk)/sqrt(2)drtrunc)/2 * 4 pi G zdisk**2 log(z/zdisk)
c  where r is spherical radius. f(r) is here taken as an exponential.
c
c  The corresponding density is:
c  f(r)*erfc*sech(z/zdisk)**2 + radial gradient terms.
c  For radii below one scale radius, we have replaced the radial exponential 
c  (which gives rise to a singular laplacian) with a quartic that joins on 
c  smoothly.
c
      r=sqrt(s*s+z*z)
c     radial truncation factors
      t=sqrt(0.5d0)*(r-outdisk)/drtrunc
      t2=t*t
      if (t.lt.-4.0) then
         eexp=0.
         eerfc=1.
      elseif (t.lt.4.0) then
         eexp=exp(-t2)/sqrt(2*pi)/drtrunc
         eerfc=0.5*erfc(t)
      else
         eexp=0
         eerfc=0
      endif
c     radial density 
c     f is radial density, f1r is f'/r, f2 is f".
      if (r.gt.0.) then
         fac1=diskconst*zdisk**2*exp(-r/rdisk)
         f=fac1*eerfc
         f1r=-fac1*(eerfc/rdisk+eexp)/r
         f2=fac1*(eerfc/rdisk2+eexp*(2/rdisk+(r-outdisk)/drtrunc**2))
      else
         fac1=diskconst*zdisk**2
         f=fac1*eerfc
         f1r=0
         f2=0
      endif
c     vertical factors
      zz=abs(z/zdisk)
      ezz=exp(-zz)
      e2zz=ezz*ezz
      tlncosh=zz+log(0.5*(1+e2zz))
      tztanh=zz*(1-e2zz)/(1+e2zz)
      tsech2=(2*ezz/(1+e2zz))**2
      appdiskdens=f2*tlncosh+2*f1r*(tztanh+tlncosh)+f*tsech2/zdisk**2

      return
      end
