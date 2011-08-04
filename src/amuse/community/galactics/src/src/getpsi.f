      function getdiskpsi(r)

      include 'commonblocks'
      common /diskestimate/ ddens(0:nmax),dpot(0:nmax),dfr(0:nmax)

      ihi=int(r/dr)+1
      if (ihi.lt.1) ihi=1
      if (ihi.gt.nr) ihi=nr
      r1=dr*(ihi-1)
      r2=dr*ihi
      t=(r-r1)/(r2-r1)
      tm1 = 1.0 - t
      if (r.eq.0.) then
         lmaxx=0
         costheta=0
      else
         costheta=z/r
         lmaxx=lmax
      endif
      getdiskpsi = (t*dpot(ihi)+ tm1*dpot(ihi-1))

      return
      end

      function gettotalpsi(rad)

      include 'commonblocks'

      gettotalpsi = 0.
      if(idiskflag.eq.1) then
         gettotalpsi = gettotalpsi + getdiskpsi(rad)
      endif
      if(ibulgeflag.eq.1) then
         gettotalpsi = gettotalpsi + sersicpot(rad)
      endif
      if(ihaloflag.eq.1) then
         gettotalpsi = gettotalpsi + gethalopsi(rad)
      endif
      return
      end

      function gethalopsi(r)

      include 'commonblocks'
      common /haloestimate/ hdens(0:nmax),hpot(0:nmax),hfr(0:nmax)

      ihi=int(r/dr)+1
      if (ihi.lt.1) ihi=1
      if (ihi.gt.nr) ihi=nr
      r1=dr*(ihi-1)
      r2=dr*ihi
      t=(r-r1)/(r2-r1)
      tm1 = 1.0 - t
      if (r.eq.0.) then
         lmaxx=0
         costheta=0
      else
         costheta=z/r
         lmaxx=lmax
      endif
      gethalopsi = (t*hpot(ihi)+ tm1*hpot(ihi-1))

      return
      end



