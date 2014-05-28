c--------------------------------
      function sigz2(r)
      parameter(nrmax=1000)
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
      common /diskpars/ sigr0, disksr, nrdisk
      common /splines/ rr(0:nrmax),fdrat(0:nrmax),drat2(0:nrmax),
     +              fszrat(0:nrmax), szrat2(0:nrmax), nrspl

c velocity dispersion is fixed so that at a height of diskz the
c density of an isothermal component has fallen by sech(1)**2=0.419974
      psizh=pot(r,zdisk)
      psi0=pot(r,0.0)
      truesigz2=(psi0-psizh)/log(0.419974)	

      call splintd(rr(0),fszrat(0),szrat2(0),nrspl+1,r,fcor)
      if(fcor.lt.0) fcor=1
      sigz2=truesigz2*fcor
!      print*,'b',fcor,psi0,psizh,r	
      return
      end
