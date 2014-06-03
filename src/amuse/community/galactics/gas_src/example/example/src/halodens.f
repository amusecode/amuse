
      function halodens(r,z)
      common /flags/ idiskflag, ibulgeflag, ihaloflag

      halodens=0
      psi=pot(r,z)
      psi0=pot(0.,0.)
      if( ihaloflag .eq. 1 ) then
         halodens = halodens + densrpsi(r,psi)
      endif
      if( ihaloflag .eq. 3 ) then
         halodens = halodens + halodfdenspsi(psi,psi0)
      endif
      return
      end

