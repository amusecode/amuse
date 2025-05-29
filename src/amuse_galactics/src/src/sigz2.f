      function sigz2(r)

      include 'commonblocks'

      psizh=pot(r,zdisk)
      psi00=pot(r,0.0)
cc
      truesigz2=(psizh-psi00)/log(0.419974)
cc
      call splintd(rr(0),fszrat(0),szrat2(0),nrspl+1,r,fcor)
      sigz2=truesigz2*fcor
      if(sigz2.lt.0.) sigz2 = 1.e-10

      return
      end
