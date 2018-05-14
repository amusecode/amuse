***
      real*8 FUNCTION mlwind(kw,lum,r,mt,mc,rl,z)
      implicit none
      integer kw
      real*8 lum,r,mt,mc,rl,z,sigma,lgdms,edfac,FOURPI,tem,ind,Xh
      real*8 ad11,ad12,ad13,ad14,ad15,ad21,ad22,ad23,ad24,ad25,ad26
      real*8 ad01,ad02,ad03,ad04,ad31,ad32,ad33,ad34
      real*8 gamma0,beta
      real*8 dml,dms,dmt,p0,x,xx,mew,neta,bwind,mxns
      real*8 lum0,kap,V1,V2,Zsun,eddington
      parameter(lum0=7.0d+04,kap=-0.5d0,V1=1.3d0,V2=2.6d0,Zsun=0.02)
      common /value1/ neta,bwind,mxns
      external eddington
*
* Calculate stellar wind mass loss.
***
* Calculate the effective temperature (tem) by the simple formula:
*     L = 4*pi*T^4*sigma
*
      sigma = 5.67d0*10.d0**(-5.d0)*(6.96d0*10.d0**10.d0)**2.d0
     &        /(3.84d0*10.d0**33.d0)
      FOURPI = 2.d0*ACOS(-1.d0)
      tem = (lum/(FOURPI*r**2.d0*sigma))**(1.d0/4.d0)
*
* Exponent of the dependence on metallicity: (Z/Zsun)^ind
* Dipendence on eddingtin factor taken from Chen (2015)
* Eddington factor is given by eddinton()
      ind = 0.85d0
      edfac = eddington(mt,lum,kw)
      if(1.d0.gt.edfac.and.edfac.gt.2.d0/3.d0)then
         ind = (2.45d0 - 2.4d0*edfac)
      elseif(edfac.gt.1.d0)then
         ind = 0.05d0
      endif
*
***
* Apply mass loss of Nieuwenhuijzen & de Jager, A&A, 1990, 231, 134,
* for massive stars over the entire HRD.
*
      dms = 0.d0
      if(lum.gt.4000.d0)then
         x = MIN(1.d0,(lum-4000.d0)/500.d0)
         dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
         dms = dms*(z/Zsun)**(ind)
      endif
*
* Stellar winds for O/B stars.
*
      if(kw.ge.0 .and. kw.le.1)then
*
* Eqs 25 and 24 in Vink et al. 2001, A&A, 369, 574.
*
          if(12500.d0.le.tem .and. tem.le.25000.d0)then
             ad11 = 2.21d0*Log10(lum/(10.d0**5.d0))
             ad12 = -1.339d0*Log10(mt/30.d0)
* V is the ration of wind velocity at inf to escape velosity.
             ad13 = -1.601d0*Log10(V1/2.d0)
             ad14 = ind*Log10(z/Zsun)
             ad15 = 1.071*Log10(tem/20000.d0)
* logarithm of the mass loss rate.
             lgdms = -6.688d0 + ad11 + ad12 + ad13 + ad14 + ad15
             dms = 10.d0**lgdms
          elseif(25000.d0.lt.tem .and. tem.le.50000.d0)then
             ad21 = 2.194d0*Log10(lum/(10.d0**5.d0))
             ad22 = -1.313d0*Log10(mt/30.d0)
* V is the ration of wind velocity at inf to escape velosity.
             ad23 = -1.226d0*Log10(V2/2.d0)
             ad24 = ind*Log10(z/Zsun)
             ad25 = 0.933d0*Log10(tem/40000.d0)
             ad26 = -10.92d0*(Log10(tem/40000.d0))**2.d0
* logarithm of the mass loss rate
             lgdms = -6.697d0 + ad21 + ad22 + ad23 + ad24 + ad25 +ad26
             dms = 10.d0**lgdms
          endif
      endif
*
***
*
      if(kw.ge.2.and.kw.le.9)then
* 'Reimers' mass loss
         dml = neta*4.0d-13*r*lum/mt
*
* Check for any tidally enhanced mass loss in binary systems (optional):
* see Tout & Eggleton, MNRAS, 1988, 231, 823.
*
         if(rl.gt.0.d0) dml = dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)

         dms = max(dms,dml)
*
* Apply mass loss of Vassiliadis & Wood, ApJ, 1993, 413, 641,
* for high pulsation periods on AGB.
*
         if(kw.eq.6)then
            p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
            p0 = 10.d0**p0
            p0 = MIN(p0,2000.d0)
            dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
            dmt = 10.d0**dmt
            dmt = 1.d0*MIN(dmt,1.36d-09*lum)
            dms = MAX(dms,dmt)
            endif
         if(kw.gt.6)then
*
***
* WR-like mass loss from Belczynski 2010 plus the dependence on 
* metallicity (Giacobbo et al. 2017).
*
            dms = (10.d0**(-13.d0))*(lum**1.5d0)*(z/Zsun)**ind
*
         endif
*
         xx = 1.0d-5*r*sqrt(lum)
         if(lum.gt.6.0d+05.and.xx.gt.1.d0)then
*
* LBV-like mass loss beyond the Humphreys-Davidson limit
* and beyond the MS (it depends on Z)
*
            dms = (z/Zsun)**ind*1.5d0*10.d0**(-4.d0)
*
         endif
      endif
*
      mlwind = dms
*
      return
      end
***
