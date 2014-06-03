      
      function diskdf3ez(ep,am,ez)
      parameter (pi=3.1415926535,nrmax=1000)
      
      psi00 = pot(0.0,0.0)
      rc=rcirc(am)
      call omekap(rc,fomega,fkappa)
      vc=rc*fomega
      ec=pot(rc,0.0)+0.5*vc**2
      if (am.lt.0.) ec = 2.0*psi00 - ec
      sr2=sigr2(rc)
      sz2=sigz2(rc)
      if (sz2.gt.0.) then
         fvert=fnamidden(rc)*exp(-ez/sz2)/sqrt(2*pi*sz2)
         diskdf3ez=fomega/(pi*fkappa)/sr2*exp(-(ep - ec)/sr2)*fvert
      else
         diskdf3ez=0
      endif
c      if (.not.(diskdf3ez.gt.0.) .and..not.(diskdf3ez.lt.0.) .and..not.(diskdf3ez.eq.0.)) then
c         write(*,*) 'DISK DF WAS NAN.'
c         write(*,*) 'EP,EZ,AM:',ep,ez,am
c         write(*,*) 'RC,EC,SR2,SZ2:',rc,ec,sr2,sz2
c         write(*,*)
c      endif
      return
      end
      
