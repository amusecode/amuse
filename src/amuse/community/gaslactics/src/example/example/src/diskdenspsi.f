
      function diskdenspsi(r,z,psi) 
      parameter (pi=3.1415926535)
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
 
c the truncation factor is derived assuming that the azimuthal velocity dist
c is gaussian with dispersion sigma_R * 2 omega/kappa, as per epicycle theory.
      if( abs(z/zdisk) .gt. 10.) goto 99
      erfarg=(r-outdisk)/1.4142136/drtrunc
      if (erfarg.gt.4.) goto 99
      if (erfarg.lt.-4) then
         trunfac=1
      else
         trunfac=0.5*erfc(erfarg)
      endif
      if (z.eq.0.) then
         con=1
      else
        psizh= pot(r,zdisk)
        psi0 = pot(r,0.)
        con = 0.419974**((psi-psi0)/(psizh-psi0))
      endif
      diskdenspsi = diskconst*exp(-r/rdisk)*con*trunfac
      return
 99   diskdenspsi = 0.0
      return
      end
