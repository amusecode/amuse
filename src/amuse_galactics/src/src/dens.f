      function dens(r,z)
      
      include 'commonblocks'

      if( idiskflag .eq. 1 ) then
         addens = appdiskdens(r,z)
         dens=totdens(r,z)
         dens = dens - addens
      else
         dens=totdens(r,z)
      endif
      return
      end
