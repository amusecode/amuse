      function diskdf5intez(vt,r,z)
      parameter (pi=3.1415926535)

      real psir0,psirz

      psir0=pot(r,0.0)
      if (z.eq.0.) then 
         psirz=psir0
      else
         psirz=pot(r,z)
         endif
      ep=0.5*(vt*vt)-psir0
      am=r*vt
      ez=-psirz+psir0
      
      diskdf5intez=diskdf3intez(ep,am,ez)
      if(diskdf5intez.lt.0.) diskdf5intez = 0.

      return
      end            


