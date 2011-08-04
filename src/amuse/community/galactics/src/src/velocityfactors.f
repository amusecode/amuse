      subroutine velocityfactors(R,z,vfacr,vfacz)

      common /blackhole/ bhmass
      real pot

      rsph = sqrt(R*R+z*z)
      call force(R,z,fr,fz,pot)
      frbh = -bhmass*R/(rsph**3.)
      fzbh = -bhmass*z/(rsph**3.)
      vfacr = sqrt(abs((fr+frbh)/fr))
      vfacz = sqrt(abs((fz+fzbh)/fz))
      if(vfacr.lt.1.or.vfacz.lt.1) then
         write(0,*) vfacr,vfacz,fr,fz
         stop
      endif
      return
      end
      
