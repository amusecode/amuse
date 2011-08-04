      function fnamidden(r)
      parameter(nrmax=1000)
      
      common /splines/ rr(0:nrmax),fdrat(0:nrmax),drat2(0:nrmax),
     +     fszrat(0:nrmax), szrat2(0:nrmax), nrspl
      
      call splintd(rr(0),fdrat(0),drat2(0),nrspl+1,r,fcor)
      fnamidden=diskdensf(r,0.0)*fcor

      return
      end
