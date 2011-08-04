      function diskdf5ez(vr,vt,vz,r,z)

      real psirz,psir0

      psir0=pot(r,0.0)
      if (z.eq.0.) then 
         psirz=psir0
      else
         psirz=pot(r,z)
      endif
c
      ep=0.5*(vr*vr+vt*vt)-psir0
c
      am=r*vt
c
      ez=0.5*vz*vz-psirz+psir0
c
      diskdf5ez=diskdf3ez(ep,am,ez)
      return
      end            


 
