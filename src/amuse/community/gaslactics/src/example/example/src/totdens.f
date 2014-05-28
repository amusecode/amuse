
      function totdens(r,z)
      common /flags/ idiskflag, ibulgeflag, ihaloflag
     
      totdens = 0
      psi=pot(r,z)
      psi0=pot(0.,0.)
!      print*,psi,diskdens(r,z,psi),bulgedenspsi(psi),densrpsi(r,psi)
      if( idiskflag .eq. 1 .or. idiskflag.eq.3) then
         totdens = totdens + diskdens(r,z,psi)
      endif
      if( ihaloflag .eq. 1 ) then
         totdens = totdens + densrpsi(r,psi)
      endif
      if( ihaloflag .eq. 3 ) then
         totdens = totdens + halodfdenspsi(psi,psi0)
      endif
      if( ibulgeflag .eq. 1 ) then
         totdens = totdens + bulgedenspsi(psi,psi0)
      endif
      if( ibulgeflag .eq. 3 ) then
         totdens = totdens + bulgedfdenspsi(psi,psi0)
      endif
      return
      end

