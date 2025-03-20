      function diskdensf(r,z)

      real psi

      psi = pot(r,z)
      diskdensf = diskdens(r,z,psi)

      return
      end
