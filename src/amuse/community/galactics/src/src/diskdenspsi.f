
      function diskdenspsi(r,z,psi) 

      include 'commonblocks'

      real psi
 
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
        psi00 = pot(r,0.)
        con = 0.419974**((psi-psi00)/(psizh-psi00))
      endif
      diskdenspsi = diskconst*exp(-r/rdisk)*con*trunfac
      return
 99   diskdenspsi = 0.0
      return
      end
