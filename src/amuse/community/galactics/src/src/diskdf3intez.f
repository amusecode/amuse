      function diskdf3intez(ep,am,ez)
      parameter (pi=3.1415926535,nrmax=1000)
      real psi00
      psi00 = pot(0.0,0.0)
      rc=rcirc(am)
      call omekap(rc,fomega,fkappa)
      vc=rc*fomega
      ec= -pot(rc,0.0)+0.5*vc*vc
      if (am.lt.0.) ec= -2*psi00-ec
      sr2=sigr2(rc)
      sz2=sigz2(rc)
      if (sz2.gt.0.) then
         fvert=fnamidden(rc)*exp(-ez/sz2)
         arg = (ep-ec)/sr2
         diskdf3intez=fomega/fkappa*sqrt(2/pi/sr2)*
     *        exp(-(arg))*fvert
      else
         diskdf3intez=0
      endif
      return
      end
