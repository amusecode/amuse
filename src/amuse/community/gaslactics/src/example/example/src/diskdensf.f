
      function diskdensf(r,z)
      psi = pot(r,z)
c      write(*,*) 'pot ',psi,r,z
      diskdensf = diskdens(r,z,psi)
      return
      end
