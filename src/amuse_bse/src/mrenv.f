***
      SUBROUTINE mrenv(kw,mass,mt,mc,lum,rad,rc,aj,tm,ltms,lbgb,lhei,
     &                 rzams,rtms,rg,menv,renv,k2e)
      implicit none
      integer kw
      real*8 mass,mt,mc,lum,rad,rc,aj,tm
      real*8 k2e,menv,menvg,menvt,menvz,renv,renvg,renvt,renvz
      real*8 A,B,C,D,E,F,x,y
      real*8 k2bgb,k2g,k2z,logm,logmt,lbgb,ltms,lhei,rg,rtms,rzams
      real*8 teff,tebgb,tetms,tau,tauenv,tautms
*
* A function to estimate the mass and radius of the convective envelope,
* as well as the gyration radius of the envelope.
* N.B. Valid only for Z=0.02!
*
* The following input is needed from HRDIAG:
*   kw = stellar type
*   mass = zero-age stellar mass
*   mt = actual mass
*   mc = core mass (not really needed, can also be done outside subroutine)
*   lum = luminosity
*   rad = radius
*   rc = core radius (not really needed...)
*   aj = age
*   tm = main-sequence lifetime
*   ltms = luminosity at TMS, lums(2)
*   lbgb = luminosity at BGB, lums(3)
*   lhei = luminosity at He ignition, lums(4)
*   rzams = radius at ZAMS
*   rtms = radius at TMS
*   rg = giant branch or Hayashi track radius, approporaite for the type. 
*        For kw=1 or 2 this is radius at BGB, and for kw=4 either GB or 
*        AGB radius at present luminosity.
*
      logm = log10(mass)
      A = MIN(0.81d0,MAX(0.68d0,0.68d0+0.4d0*logm))
      C = MAX(-2.5d0,MIN(-1.5d0,-2.5d0+5.d0*logm))
      D = -0.1d0
      E = 0.025d0
*
* Zero-age and BGB values of k^2.
*
      k2z = MIN(0.21d0,MAX(0.09d0-0.27d0*logm,0.037d0+0.033d0*logm))
      if(logm.gt.1.3d0) k2z = k2z - 0.055d0*(logm-1.3d0)**2
      k2bgb = MIN(0.15d0,MIN(0.147d0+0.03d0*logm,0.162d0-0.04d0*logm))
*
      if(kw.ge.3.and.kw.le.6)then
*
* Envelope k^2 for giant-like stars; this will be modified for non-giant
* CHeB stars or small envelope mass below.
* Formula is fairly accurate for both FGB and AGB stars if M <= 10, and
* gives reasonable values for higher masses. Mass dependence is on actual
* rather than ZA mass, expected to work for mass-losing stars (but not
* tested!). The slightly complex appearance is to insure continuity at 
* the BGB, which depends on the ZA mass.
*
         logmt = log10(mt)
         F = 0.208d0 + 0.125d0*logmt - 0.035d0*logmt**2
         B = 1.0d+04*mt**(3.d0/2.d0)/(1.d0+0.1d0*mt**(3.d0/2.d0))
         x = ((lum-lbgb)/B)**2
         y = (F - 0.033d0*log10(lbgb))/k2bgb - 1.d0
         k2g = (F - 0.033d0*log10(lum) + 0.4d0*x)/(1.d0+y*(lbgb/lum)+x)
      elseif(kw.eq.9)then
*
* Rough fit for for HeGB stars...
*
         B = 3.0d+04*mt**(3.d0/2.d0)
         x = (MAX(0.d0,lum/B-0.5d0))**2
         k2g = (k2bgb + 0.4d0*x)/(1.d0 + 0.4d0*x)
      else
         k2g = k2bgb
      endif
*
      if(kw.le.2)then
         menvg = 0.5d0
         renvg = 0.65d0
      elseif(kw.eq.3.and.lum.lt.3.d0*lbgb)then
*
* FGB stars still close to the BGB do not yet have a fully developed CE.
*
         x = MIN(3.d0,lhei/lbgb)
         tau = MAX(0.d0,MIN(1.d0,(x-lum/lbgb)/(x-1.d0)))
         menvg = 1.d0 - 0.5d0*tau**2
         renvg = 1.d0 - 0.35d0*tau**2
      else
         menvg = 1.d0
         renvg = 1.d0
      endif
*
      if(rad.lt.rg)then
*
* Stars not on the Hayashi track: MS and HG stars, non-giant CHeB stars,
* HeMS and HeHG stars, as well as giants with very small envelope mass.
*
         
         if(kw.le.6)then
*
* Envelope k^2 fitted for MS and HG stars.
* Again, pretty accurate for M <= 10 but less so for larger masses.
* [Note that this represents the whole star on the MS, so there is a 
* discontinuity in stellar k^2 between MS and HG - okay for stars with a 
* MS hook but low-mass stars should preferably be continous...]
*
* For other types of star not on the Hayashi track we use the same fit as 
* for HG stars, this is not very accurate but has the correct qualitative 
* behaviour. For CheB stars this is an overestimate because they appear
* to have a more centrally concentrated envelope than HG stars.
*
            k2e = (k2z-E)*(rad/rzams)**C + E*(rad/rzams)**D
         elseif(kw.eq.7)then
* Rough fit for naked He MS stars.
            tau = aj/tm
            k2e = 0.08d0 - 0.03d0*tau
         elseif(kw.le.9)then
* Rough fit for HeHG stars.
            k2e = 0.08d0*rzams/rad
         endif
*
* tauenv measures proximity to the Hayashi track in terms of Teff.
* If tauenv>0 then an appreciable convective envelope is present, and
* k^2 needs to be modified.
*
         if(kw.le.2)then
            teff = sqrt(sqrt(lum)/rad)
            tebgb = sqrt(sqrt(lbgb)/rg)
            tauenv = MAX(0.d0,MIN(1.d0,(tebgb/teff-A)/(1.d0-A)))
         else
            tauenv = MAX(0.d0,MIN(1.d0,(sqrt(rad/rg)-A)/(1.d0-A)))
         endif
*
         if(tauenv.gt.0.d0)then
            menv = menvg*tauenv**5
            renv = renvg*tauenv**(5.d0/4.d0)
            if(kw.le.1)then
* Zero-age values for CE mass and radius.
               x = MAX(0.d0,MIN(1.d0,(0.1d0-logm)/0.55d0))
               menvz = 0.18d0*x + 0.82d0*x**5
               renvz = 0.4d0*x**(1.d0/4.d0) + 0.6d0*x**10
               y = 2.d0 + 8.d0*x
* Values for CE mass and radius at start of the HG.
               tetms = sqrt(sqrt(ltms)/rtms)
               tautms = MAX(0.d0,MIN(1.d0,(tebgb/tetms-A)/(1.d0-A)))
               menvt = menvg*tautms**5
               renvt = renvg*tautms**(5.d0/4.d0)
* Modified expressions during MS evolution.
               tau = aj/tm
               if(tautms.gt.0.d0)then
                  menv = menvz + tau**y*menv*(menvt - menvz)/menvt
                  renv = renvz + tau**y*renv*(renvt - renvz)/renvt
               else
                  menv = 0.d0
                  renv = 0.d0
               endif
               k2e = k2e + tau**y*tauenv**3*(k2g - k2e)
            else
               k2e = k2e + tauenv**3*(k2g - k2e)
            endif
         else
            menv = 0.d0
            renv = 0.d0
         endif
      else
*
* All other stars should be true giants.
*
         menv = menvg
         renv = renvg
         k2e = k2g
      endif
*
      menv = menv*(mt - mc)
      renv = renv*(rad - rc)
      menv = MAX(menv,1.0d-10)
      renv = MAX(renv,1.0d-10)
*
      return
      end
***
