      function sigr2(r)

      common /diskpars/ sigr0, disksr, nrdisk

      sigr2=sigr0*sigr0*exp(-r/disksr)

      return
      end
