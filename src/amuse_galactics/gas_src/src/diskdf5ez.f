c--------------------
      function diskdf5ez(vr,vt,vz,r,z)
      parameter (pi=3.1415926535)
      parameter(nrmax=1000)

      psir0=pot(r,0.0)
      if (z.eq.0.) then 
         psirz=psir0
      else
         psirz=pot(r,z)
         endif
      ep=0.5*(vr*vr+vt*vt)+psir0
      am=r*vt
      ez=0.5*vz*vz+psirz-psir0
      diskdf5ez=diskdf3ez(ep,am,ez)
      return
      end            
