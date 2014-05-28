      function pot(s,z)
      parameter (pi=3.1415926535)
      common /potconstants/ apot(20,0:20000), frad(20,0:20000), 
     +     dr, nr, lmax, potcor
      common /legendre/ plcon(0:40)
      common /flags/ idiskflag, ibulgeflag, ihaloflag
c     Returns the potential, being the sum of the spherical harmonics
c     with coefficients in apot, and the potential appdiskpot which
c     approximates the high-frequency components of the disk.
c     + external potentials      
      p=0
      
      r=sqrt(s*s+z*z)
      ihi=int(r/dr)+1
      if (ihi.lt.1) ihi=1
      if (ihi.gt.nr) ihi=nr
      r1=dr*(ihi-1)
      r2=dr*ihi
      t=(r-r1)/(r2-r1)
      tm1 = 1.0 - t
      if (r.eq.0.) then
         lmaxx=0
         costheta=0
      else
         costheta=z/r
         lmaxx=lmax
      endif
      
      do l=lmaxx,0,-2
         p=p+plgndr1(l,costheta)*plcon(l)*
     +        (t*apot(l/2+1,ihi)+ tm1*apot(l/2+1,ihi-1))
      enddo
      if( idiskflag .eq. 1 .or. idiskflag .eq. 3)
     & p = p + appdiskpot(s,z)
      
      if( idiskflag .eq. 2 .or. idiskflag .eq. 3)
     & p= p + thickdiskpot(s,z)      

      if(ihaloflag.eq.2)
     & p= p + fixedpot(r)
      
      pot=p
      return
      end

      function extpot(s,z)
      parameter (pi=3.1415926535)
      common /flags/ idiskflag, ibulgeflag, ihaloflag
c     external potentials      
      p=0
      r=sqrt(s*s+z*z) 
           
      if( idiskflag .eq. 1 .or. idiskflag .eq. 3)
     & p = p + appdiskpot(s,z)

      if( idiskflag .eq. 2 .or. idiskflag .eq. 3)
     & p= p + thickdiskpot(s,z)      

      if(ihaloflag.eq.2)
     & p= p + fixedpot(r)
      
      extpot=p
      return
      end



      function cfharm(s,z,ff)
      parameter (pi=3.1415926535)
      common /potconstants/ apot(20,0:20000), frad(20,0:20000), 
     +     dr, nr, lmax, potcor
      common /legendre/ plcon(0:40)
      common /flags/ idiskflag, ibulgeflag, ihaloflag
c     Returns the potential, being the sum of the spherical harmonics
c     with coefficients in apot, and the potential appdiskpot which
c     approximates the high-frequency components of the disk.
c     + external potentials
      real ff(20,0:20000)      
      p=0
      r=sqrt(s*s+z*z)
      ihi=int(r/dr)+1
      if (ihi.lt.1) ihi=1
      if (ihi.gt.nr) ihi=nr
      r1=dr*(ihi-1)
      r2=dr*ihi
      t=(r-r1)/(r2-r1)
      tm1 = 1.0 - t
      if (r.eq.0.) then
         lmaxx=0
         costheta=0
      else
         costheta=z/r
         lmaxx=lmax
      endif
      
      do l=lmaxx,0,-2
         p=p+plgndr1(l,costheta)*plcon(l)*
     +        (t*ff(l/2+1,ihi)+ tm1*ff(l/2+1,ihi-1))
      enddo
            
      cfharm=p
      return
      end
