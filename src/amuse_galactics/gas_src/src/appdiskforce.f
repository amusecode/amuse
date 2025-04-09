
      subroutine appdiskforce(s,z,fsad,fzad)
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea

      parameter (pi=3.1415926535)

c      fsad = 0
c      fzad = 0
c      return
      r=sqrt(s*s + z*z)
      t=sqrt(0.5d0)*(r-outdisk)/drtrunc
      if (t.lt.-4.0) then
          eerfc=1.
          derfc=0.0
      elseif(t.gt.4.0) then
          eerfc = 0.0
          derfc = 0.0
      else
          eerfc=0.5*erfc(t)
          derfc = exp(-t**2)/(sqrt(2.0*pi)*drtrunc)
      endif

      r1 = r/rdisk
      if (r1.eq.0.) then
         fprime=0
         f=0
      else
          texp = diskconst*zdisk**2*exp(-r1)
          f = texp*eerfc
          fprime = -texp*(derfc + eerfc/rdisk)/r
      endif

      zz = abs(z/zdisk)
      if( zz .lt. 10 ) then
          ezz = exp(-2.0*zz)
      else
          e2zz = 0.0
      endif
      tlncoshz = zz + log(0.5*(1.0 + e2zz))

      fsad = -4.0*pi*fprime*s*tlncoshz
      fzad = -4.0*pi*(fprime*z*tlncoshz + f/zdisk*tanh(z/zdisk))

      return
      end
