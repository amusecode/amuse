***
      real*8 FUNCTION mlwind(kw,lum,r,mt,mc,rl,z)
      implicit none
      integer kw
      real*8 lum,r,mt,mc,rl,z
      real*8 dml,dms,dmt,p0,x,mew,lum0,kap,neta,bwind,hewind,mxns
      parameter(lum0=7.0d+04,kap=-0.5d0)
      common /value1/ neta,bwind,hewind,mxns
*
* Calculate stellar wind mass loss.
*
* Apply mass loss of Nieuwenhuijzen & de Jager, A&A, 1990, 231, 134,
* for massive stars over the entire HRD.
      dms = 0.d0
      if(lum.gt.4000.d0)then
         x = MIN(1.d0,(lum-4000.d0)/500.d0)
         dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
         dms = dms*(z/0.02d0)**(1.d0/2.d0)
      endif
      if(kw.ge.2.and.kw.le.9)then
* 'Reimers' mass loss
         dml = neta*4.0d-13*r*lum/mt
         if(rl.gt.0.d0) dml = dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
* Apply mass loss of Vassiliadis & Wood, ApJ, 1993, 413, 641, 
* for high pulsation periods on AGB.
         if(kw.eq.5.or.kw.eq.6)then
            p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
            p0 = 10.d0**p0
            p0 = MIN(p0,2000.d0)
            dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
            dmt = 10.d0**dmt
            dmt = 1.d0*MIN(dmt,1.36d-09*lum)
            dml = MAX(dml,dmt)
         endif
         if(kw.gt.6)then
            dms = MAX(dml,1.0d-13*hewind*lum**(3.d0/2.d0))
         else
            dms = MAX(dml,dms)
            mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
* reduced WR-like mass loss for small H-envelope mass
            if(mew.lt.1.d0)then
               dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
               dms = MAX(dml,dms)
            end if
* LBV-like mass loss beyond the Humphreys-Davidson limit.
            x = 1.0d-5*r*sqrt(lum)
            if(lum.gt.6.0d+05.and.x.gt.1.d0)then
               dml = 0.1d0*(x-1.d0)**3*(lum/6.0d+05-1.d0)
               dms = dms + dml
            endif
         endif
      endif
*
      mlwind = dms
*
      return
      end
***
