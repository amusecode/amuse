***
      SUBROUTINE hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r,lum,kw,mc,rc,menv,renv,k2)
*
*
*       H-R diagram for population I stars.
*       -----------------------------------
*
*       Computes the new mass, luminosity, radius & stellar type.
*       Input (MASS, AJ, TM, TN, LUMS & TSCLS) supplied by routine STAR.
*       Ref: P.P. Eggleton, M.J. Fitchett & C.A. Tout (1989) Ap.J. 347, 998.
*
*       Revised 27th March 1995 by C. A. Tout;
*       24th October 1995 to include metallicity;
*       14th November 1996 to include naked helium stars;
*       28th February 1997 to allow accretion induced supernovae.
*
*       Revised 5th April 1997 by J. R. Hurley
*       to include from Z=0.0001 to Z=0.03, convective overshooting,
*       MS hook and more elaborate CHeB.
*
      implicit none
*
      integer kw,kwp
      integer wdflag,nsflag
      parameter(wdflag=1,nsflag=1)
      common /fall/fallback
*
      real*8 mass,aj,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
      real*8 r,lum,mc,rc,menv,renv,k2,fallback
      real*8 mch,mlp,tiny
      parameter(mch=1.44d0,mlp=12.d0,tiny=1.0d-14)
      real*8 mxns,mxns0,mxns1
      parameter(mxns0=1.8d0,mxns1=3.d0)
      real*8 mass0,mt0,mtc
* 
      real*8 thook,thg,tbagb,tau,tloop,taul,tauh,tau1,tau2,dtau,texp
      real*8 lx,ly,dell,alpha,beta,neta
      real*8 rx,ry,delr,rzams,rtms,gamma,rmin,taumin,rg
      parameter(taumin=5.0d-08)
      real*8 mcmax,mcx,mcy,mcbagb,lambda
      real*8 am,xx,fac,rdgen,mew,lum0,kap,zeta,ahe,aco
      parameter(lum0=7.0d+04,kap=-0.5d0,ahe=4.d0,aco=16.d0)
*
      real*8 thookf,tblf
      real*8 lalphf,lbetaf,lnetaf,lhookf,lgbtf,lmcgbf,lzhef,lpertf
      real*8 rzamsf,rtmsf,ralphf,rbetaf,rgammf,rhookf
      real*8 rgbf,rminf,ragbf,rzahbf,rzhef,rhehgf,rhegbf,rpertf
      real*8 mctmsf,mcgbtf,mcgbf,mcheif,mcagbf,lzahbf
      external thookf,tblf
      external lalphf,lbetaf,lnetaf,lhookf,lgbtf,lmcgbf,lzhef,lpertf
      external rzamsf,rtmsf,ralphf,rbetaf,rgammf,rhookf
      external rgbf,rminf,ragbf,rzahbf,rzhef,rhehgf,rhegbf,rpertf
      external mctmsf,mcgbtf,mcgbf,mcheif,mcagbf,lzahbf
*
*
*       ---------------------------------------------------------------------
*       MASS    Stellar mass in solar units (input: old; output: new value).
*       AJ      Current age in Myr.
*       MT      Current mass in solar units (used for R).
*       TM      Main sequence time.
*       TN      Nuclear burning time.
*       TSCLS   Time scale for different stages.
*       LUMS    Characteristic luminosity.
*       GB      Giant Branch parameters
*       ZPARS   Parameters for distinguishing various mass intervals.
*       R       Stellar radius in solar units.
*       TE      Effective temperature (suppressed).
*       KW      Classification type (0 - 15).
*       MC      Core mass.
*       fallback Belczynski's Prescription for fallback on black holes
*       ---------------------------------------------------------------------
*
      fallback = 0.d0
* Set maximum NS mass depending on which NS mass prescription is used. 
      mxns = mxns0
      if(nsflag.eq.1) mxns = mxns1
*
      mass0 = mass
      if(mass0.gt.100.d0) mass = 100.d0
      mt0 = mt
      if(mt0.gt.100.d0) mt = 100.d0
*
* Make evolutionary changes to stars that have not reached KW > 5.
*
      if(kw.gt.6) goto 90
*
      tbagb = tscls(2) + tscls(3)
      thg = tscls(1) - tm
*
      rzams = rzamsf(mass)
      rtms = rtmsf(mass)
*
      if(aj.lt.tscls(1))then
*
*        Either on MS or HG
*
         rg = rgbf(mt,lums(3))
*
         if(aj.lt.tm)then
*
*           Main sequence star.
*
            mc = 0.d0
            tau = aj/tm
            thook = thookf(mass)*tscls(1)
            zeta = 0.01d0
            tau1 = MIN(1.d0,aj/thook)
            tau2 = MAX(0.d0,
     &             MIN(1.d0,(aj-(1.d0-zeta)*thook)/(zeta*thook)))
*
            dell = lhookf(mass,zpars(1))
            dtau = tau1**2 - tau2**2
            alpha = lalphf(mass)
            beta = lbetaf(mass)
            neta = lnetaf(mass)
            lx = LOG10(lums(2)/lums(1))
            if(tau.gt.taumin)then
               xx = alpha*tau + beta*tau**neta +
     &              (lx - alpha - beta)*tau**2 - dell*dtau
            else
               xx = alpha*tau + (lx - alpha)*tau**2 - dell*dtau
            endif
            lum = lums(1)*10.d0**xx
*
            delr = rhookf(mass,zpars(1))
            dtau = tau1**3 - tau2**3
            alpha = ralphf(mass)
            beta = rbetaf(mass)
            gamma = rgammf(mass)
            rx = LOG10(rtms/rzams)
* Note that the use of taumin is a slightly pedantic attempt to
* avoid floating point underflow. It IS overkill!
            if(tau.gt.taumin)then
               xx = alpha*tau + beta*tau**10 + gamma*tau**40 +
     &              (rx - alpha - beta - gamma)*tau**3 - delr*dtau
            else
               xx = alpha*tau + (rx - alpha)*tau**3 - delr*dtau
            endif
            r = rzams*10.d0**xx
*
            if(mass.lt.(zpars(1)-0.3d0))then
               kw = 0
* This following is given by Chris for low mass MS stars which will be 
* substantially degenerate. We need the Hydrogen abundance, X, which we 
* calculate from Z assuming that the helium abundance, Y, is calculated
* according to Y = 0.24 + 2*Z
               rdgen = 0.0258d0*((1.d0+zpars(11))**(5.d0/3.d0))*
     &                          (mass**(-1.d0/3.d0))
               r = MAX(rdgen,r)
            else
               kw = 1
            endif
*
         else 
*
*           Star is on the HG
*
            mcx = mc
            if(mass.le.zpars(2))then
               mc = mcgbf(lums(3),GB,lums(6))
            elseif(mass.le.zpars(3))then
               mc = mcheif(mass,zpars(2),zpars(9))
            else
               mc = mcheif(mass,zpars(2),zpars(10))
            endif
            neta = mctmsf(mass)
            tau = (aj - tm)/thg
            mc = ((1.d0 - tau)*neta + tau)*mc
            mc = MAX(mc,mcx)
*
* Test whether core mass has reached total mass.
*
            if(mc.ge.mt)then
               aj = 0.d0
               if(mass.gt.zpars(2))then
*
* Zero-age helium star
*
                  mc = 0.d0
                  mass = mt
                  kw = 7
                  CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               else
*
* Zero-age helium white dwarf.
*
                  mc = mt
                  mass = mt
                  kw = 10
               endif
            else
               lum = lums(2)*(lums(3)/lums(2))**tau
               if(mass.le.zpars(3))then
                  rx = rg
               else
* He-ignition and end of HG occur at Rmin
                  rmin = rminf(mass)
                  ry = ragbf(mt,lums(4),zpars(2))
                  rx = MIN(rmin,ry)
                  if(mass.le.mlp)then
                     texp = log(mass/mlp)/log(zpars(3)/mlp)
                     rx = rg
                     rx = rmin*(rx/rmin)**texp
                  endif
                  tau2 = tblf(mass,zpars(2),zpars(3))
                  if(tau2.lt.tiny) rx = ry
               endif
               r = rtms*(rx/rtms)**tau
               kw = 2
            endif
*
         endif
*
* Now the GB, CHeB and AGB evolution.
*
      elseif(aj.lt.tscls(2))then
*
*        Red Giant.
*
         kw = 3
         lum = lgbtf(aj,GB(1),GB,tscls(4),tscls(5),tscls(6))
         if(mass.le.zpars(2))then
* Star has a degenerate He core which grows on the GB
            mc = mcgbf(lum,GB,lums(6))
         else
* Star has a non-degenerate He core which may grow, but
* only slightly, on the GB
            tau = (aj - tscls(1))/(tscls(2) - tscls(1))
            mcx = mcheif(mass,zpars(2),zpars(9))
            mcy = mcheif(mass,zpars(2),zpars(10))
            mc = mcx + (mcy - mcx)*tau
         endif
         r = rgbf(mt,lum)
         rg = r
         if(mc.ge.mt)then
            aj = 0.d0
            if(mass.gt.zpars(2))then
*
* Zero-age helium star
*
               mc = 0.d0
               mass = mt
               kw = 7
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            else
*
* Zero-age helium white dwarf.
*
               mc = mt
               mass = mt
               kw = 10
            endif
         endif
*
      elseif(aj.lt.tbagb)then
*
*       Core helium burning star.
*
         if(kw.eq.3.and.mass.le.zpars(2))then
            mass = mt
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            aj = tscls(2)
         endif
         if(mass.le.zpars(2))then
            mcx = mcgbf(lums(4),GB,lums(6))
         else
            mcx = mcheif(mass,zpars(2),zpars(10))
         end if
         tau = (aj - tscls(2))/tscls(3)
         mc = mcx + (mcagbf(mass) - mcx)*tau
*
         if(mass.le.zpars(2))then
            lx = lums(5)
            ly = lums(7)
            rx = rzahbf(mt,mc,zpars(2))
            rg = rgbf(mt,lx)
            rmin = rg*zpars(13)**(mass/zpars(2))
            texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
            ry = ragbf(mt,ly,zpars(2))
            if(rmin.lt.rx)then
               taul = (log(rx/rmin))**(1.d0/3.d0)
            else
               rmin = rx
               taul = 0.d0
            endif
            tauh = (log(ry/rmin))**(1.d0/3.d0)
            tau2 = taul*(tau - 1.d0) + tauh*tau
            r = rmin*exp(abs(tau2)**3)
            rg = rg + tau*(ry - rg)
            lum = lx*(ly/lx)**(tau**texp)
         elseif(mass.gt.zpars(3))then
*
* For HM stars He-ignition takes place at Rmin in the HG, and CHeB
* consists of a blue phase (before tloop) and a RG phase (after tloop).
*
            tau2 = tblf(mass,zpars(2),zpars(3))
            tloop = tscls(2) + tau2*tscls(3)
            rmin = rminf(mass)
            rg = rgbf(mt,lums(4))
            rx = ragbf(mt,lums(4),zpars(2))
            rmin = MIN(rmin, rx)
            if(mass.le.mlp) then
               texp = log(mass/mlp)/log(zpars(3)/mlp)
               rx = rg
               rx = rmin*(rx/rmin)**texp
            else
               rx = rmin
            end if
            texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
            lum = lums(4)*(lums(7)/lums(4))**(tau**texp)
            if(aj.lt.tloop)then
               ly = lums(4)*(lums(7)/lums(4))**(tau2**texp)
               ry = ragbf(mt,ly,zpars(2))
               taul = 0.d0
               if(ABS(rmin-rx).gt.tiny)then
                  taul = (log(rx/rmin))**(1.d0/3.d0)
               endif
               tauh = 0.d0
               if(ry.gt.rmin) tauh = (log(ry/rmin))**(1.d0/3.d0)
               tau = (aj - tscls(2))/(tau2*tscls(3))
               tau2 = taul*(tau - 1.d0) + tauh*tau
               r = rmin*exp(abs(tau2)**3)
               rg = rg + tau*(ry - rg)
            else
               r = ragbf(mt,lum,zpars(2))
               rg = r
            end if
         else
*
* For IM stars CHeB consists of a RG phase (before tloop) and a blue
* loop (after tloop).
*
            tau2 = 1.d0 - tblf(mass,zpars(2),zpars(3))
            tloop = tscls(2) + tau2*tscls(3)
            if(aj.lt.tloop)then
               tau = (tloop - aj)/(tau2*tscls(3))
               lum = lums(5)*(lums(4)/lums(5))**(tau**3)
               r = rgbf(mt,lum)
               rg = r
            else
               lx = lums(5)
               ly = lums(7)
               rx = rgbf(mt,lx)
               rmin = rminf(mt)
               texp = MIN(MAX(0.4d0,rmin/rx),2.5d0)
               ry = ragbf(mt,ly,zpars(2))
               if(rmin.lt.rx)then
                  taul = (log(rx/rmin))**(1.d0/3.d0)
               else
                  rmin = rx
                  taul = 0.d0
               endif
               tauh = (log(ry/rmin))**(1.d0/3.d0)
               tau = (aj - tloop)/(tscls(3) - (tloop - tscls(2)))
               tau2 = taul*(tau - 1.d0) + tauh*tau
               r = rmin*exp(abs(tau2)**3)
               rg = rx + tau*(ry - rx)
               lum = lx*(ly/lx)**(tau**texp)
            endif
         endif
* 
* Test whether core mass exceeds total mass.
*
         if(mc.ge.mt)then
*
* Evolved MS naked helium star.
*
            kw = 7
            xx = (aj - tscls(2))/tscls(3)
            mass = mt
            CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
            aj = xx*tm
         else
            kw = 4
         endif
*
      else
*
*        Asymptotic Red Giant.
*
* On the AGB the He core mass remains constant until at Ltp it
* is caught by the C core mass and they grow together.
*
         mcbagb = mcagbf(mass)
         mcx = mcgbtf(tbagb,GB(8),GB,tscls(7),tscls(8),tscls(9))
         mcmax = MAX(MAX(mch,0.773d0*mcbagb-0.35d0),1.05d0*mcx)
*
         if(aj.lt.tscls(13))then
            mcx = mcgbtf(aj,GB(8),GB,tscls(7),tscls(8),tscls(9))
            mc = mcbagb
            lum = lmcgbf(mcx,GB)
            if(mt.le.mc)then
*
* Evolved naked helium star as the envelope is lost but the
* star has not completed its interior burning. The star becomes
* a post-HeMS star.
*
               kw = 9
               mt = mc
               mass = mt
               mc = mcx
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               if(mc.le.GB(7))then
                  aj = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                            (mc**(1.d0-GB(5)))
               else
                  aj = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                            (mc**(1.d0-GB(6)))
               endif
               aj = MAX(aj,tm)
               goto 90
            else
               kw = 5
            endif
         else
            kw = 6
            mc = mcgbtf(aj,GB(2),GB,tscls(10),tscls(11),tscls(12))
            lum = lmcgbf(mc,GB)
*
* Approximate 3rd Dredge-up on AGB by limiting Mc.
*
            lambda = MIN(0.9d0,0.3d0+0.001d0*mass**5)
            tau = tscls(13)
            mcx = mcgbtf(tau,GB(2),GB,tscls(10),tscls(11),tscls(12))
            mcy = mc
            mc = mc - lambda*(mcy-mcx)
            mcx = mc
            mcmax = MIN(mt,mcmax)   
         endif
         r = ragbf(mt,lum,zpars(2))
         rg = r
*
* Mc,x represents the C core mass and we now test whether it
* exceeds either the total mass or the maximum allowed core mass.
*
         if(mcmax-mcx.lt.tiny)then
            aj = 0.d0
            mc = mcmax
            if(mc.lt.mch)then
               mt = mc
               if(mcbagb.lt.1.6d0)then
*     
* Zero-age Carbon/Oxygen White Dwarf
*
                  kw = 11
               else
*     
* Zero-age Oxygen/Neon White Dwarf
*
                  kw = 12
               endif
               mass = mt
            else
               if(mcbagb.lt.1.6d0)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
                  kw = 15
                  aj = 0.d0
                  mt = 0.d0
                  lum = 1.0d-10
                  r = 1.0d-10
               else
                  if(nsflag.eq.0)then
                     mt = 1.17d0 + 0.09d0*mc
                  elseif(nsflag.ge.1)then
*
* Use NS/BH mass given by Belczynski et al. 2002, ApJ, 572, 407.
*
                     if(mc.lt.2.5d0)then
                        mcx = 0.161767d0*mc + 1.067055d0
                     else
                        mcx = 0.314154d0*mc + 0.686088d0
                     endif
                     if(mc.le.5.d0)then
                        mt = mcx
                        fallback = 0.d0
                     elseif(mc.lt.7.6d0)then
                        mt = mcx + (mc - 5.d0)*(mt - mcx)/2.6d0
                        fallback = (mc - 5.d0)/2.6d0
                     endif
                     if(mc.gt.7.60) fallback = 1.d0
                  endif
                  mc = mt
                  if(mt.le.mxns)then
*
* Zero-age Neutron star
*
                     kw = 13
                  else
*
* Zero-age Black hole
*
                     kw = 14
*       print*,' HRDIAG BH mc, mt, mxns, fallback ',mc,mt,mxns,fallback
                  endif  
               endif
            endif
         endif
*
      endif
*
 90   continue
*
      if(kw.ge.7.and.kw.le.9)then
*
* Naked Helium Star
*
         rzams = rzhef(mt)
         rx = rzams
         if(aj.lt.tm)then
*
* Main Sequence
*
            kw = 7
            tau = aj/tm
            am = MAX(0.d0,0.85d0-0.08d0*mass)
            lum = lums(1)*(1.d0+0.45d0*tau+am*tau**2)
            am = MAX(0.d0,0.4d0-0.22d0*LOG10(mt))
            r = rx*(1.d0+am*(tau-tau**6))
            rg = rx
* Star has no core mass and hence no memory of its past
* which is why we subject mass and mt to mass loss for
* this phase.
            mc = 0.d0
            if(mt.lt.zpars(10)) kw = 10
         else
*
* Helium Shell Burning
*
            kw = 8
            lum = lgbtf(aj,GB(8),GB,tscls(4),tscls(5),tscls(6))
            r = rhehgf(mt,lum,rx,lums(2))
            rg = rhegbf(lum)
            if(r.ge.rg)then
               kw = 9
               r = rg
            endif
            mc = mcgbf(lum,GB,lums(6))
            mtc = MIN(mt,1.45d0*mt-0.31d0)
            mcmax = MIN(mtc,MAX(mch,0.773d0*mass-0.35d0))
            if(mcmax-mc.lt.tiny)then
               aj = 0.d0
               mc = mcmax
               if(mc.lt.mch)then
                  if(mass.lt.1.6d0)then
*     
* Zero-age Carbon/Oxygen White Dwarf
*
                     mt = MAX(mc,(mc+0.31d0)/1.45d0)
                     kw = 11
                  else
*     
* Zero-age Oxygen/Neon White Dwarf
*
                     mt = mc
                     kw = 12
                  endif
                  mass = mt
               else
                  if(mass.lt.1.6d0)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
                     kw = 15
                     aj = 0.d0
                     mt = 0.d0
                     lum = 1.0d-10
                     r = 1.0d-10
                  else
                     if(nsflag.eq.0)then
                        mt = 1.17d0 + 0.09d0*mc
                     elseif(nsflag.ge.1)then
                        if(mc.lt.2.5d0)then
                           mcx = 0.161767d0*mc + 1.067055d0
                        else
                           mcx = 0.314154d0*mc + 0.686088d0
                        endif
                        if(mc.le.5.d0)then
                           mt = mcx
                           fallback = 0.d0
                        elseif(mc.lt.7.6d0)then
                           mt = mcx + (mc - 5.d0)*(mt - mcx)/2.6d0
                           fallback = (mc - 5.d0)/2.6d0
                        endif
                        if(mc.gt.7.60) fallback = 1.d0
                     endif
                     mc = mt
                     if(mt.le.mxns)then
*
* Zero-age Neutron star
*
                        kw = 13
                     else
*
* Zero-age Black hole
*
                        kw = 14
*                     print*,' HRDIAG BH T ',mc,mt,mxns,fallback
                     endif
                  endif  
               endif
            endif
         endif
      endif
*
      if(kw.ge.10.and.kw.le.12)then
*
*        White dwarf.
*
         mc = mt
         if(mc.ge.mch.and.kw.eq.12)then
*
* Accretion induced supernova with no remnant
* unless WD is ONe.
*
            kw = 13
            mt = 1.3d0
*           mt = 1.17d0 + 0.09d0*mc
         elseif(mc.ge.(mch-0.06d0).and.kw.le.11)then
            kw = 15
            aj = 0.d0
            mt = 0.d0
            lum = 1.0d-10
            r = 1.0d-10
         else
*
            if(kw.eq.10)then
               xx = ahe
            else
               xx = aco
            endif
*
            if(wdflag.eq.0)then
*
* Mestel cooling
*
               lum = 635.d0*mt*zpars(14)/(xx*(aj+0.1d0))**1.4d0
*
            elseif(wdflag.ge.1)then
*
* modified-Mestel cooling
*
               if(aj.lt.9000.0)then
                  lum = 300.d0*mt*zpars(14)/(xx*(aj+0.1d0))**1.18d0
               else
                  fac = (9000.1d0*xx)**5.3d0
                  lum = 300.d0*fac*mt*zpars(14)/(xx*(aj+0.1d0))**6.48d0
               endif
*
            endif
*
            if(mt.lt.0.000005d0)then
               r = 0.009d0
            elseif(mt.lt.0.0005d0)then
               r = 0.09d0
            else
               r = 0.0115d0*SQRT(MAX(1.48204d-06,(mch/mt)**(2.d0/3.d0)
     &                                      - (mt/mch)**(2.d0/3.d0)))
               r = MIN(0.1d0,r)
            endif
*
         endif
      endif
*
      if(kw.eq.13)then
*
*        Neutron Star.
*
         mc = mt
         if(mc.gt.mxns)then
*
* Accretion induced Black Hole?
*
            kw = 14
            aj = 0.d0
         else
            lum = 0.02d0*(mt**0.67d0)/(MAX(aj,0.1d0))**2
            r = 1.4d-05
         endif
      endif
*
      if(kw.eq.14)then
*
*        Black hole
*
         mc = mt
         lum = 1.0d-10
         r = 4.24d-06*mt
      endif
*
* Calculate the core radius and the luminosity and radius of the
* remnant that the star will become.
*
      tau = 0.d0
      if(kw.le.1.or.kw.eq.7)then
         rc = 0.d0
      elseif(kw.le.3)then
         if(mass.gt.zpars(2))then
            lx = lzhef(mc)
            rx = rzhef(mc)
            rc = rx
         else
            if(wdflag.eq.0)then
               lx = 635.d0*mc*zpars(14)/((ahe*0.1d0)**1.4d0)
            elseif(wdflag.ge.1)then
               lx = 300.d0*mc*zpars(14)/((ahe*0.1d0)**1.18d0)
            endif
            rx = 0.0115d0*SQRT(MAX(1.48204d-06,
     &           (mch/mc)**(2.d0/3.d0)-(mc/mch)**(2.d0/3.d0)))
            rc = 5.d0*rx
         endif
      elseif(kw.eq.4)then
         tau = (aj - tscls(2))/tscls(3)
         kwp = 7
         CALL star(kwp,mc,mc,tm,tn,tscls,lums,GB,zpars)
         am = MAX(0.d0,0.85d0-0.08d0*mc)
         lx = lums(1)*(1.d0+0.45d0*tau+am*tau**2)
         rx = rzhef(mc)
         am = MAX(0.d0,0.4d0-0.22d0*LOG10(mc))
         rx = rx*(1.d0+am*(tau-tau**6))
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         rc = rx
      elseif(kw.eq.5)then
         kwp = 9
         if(tn.gt.tbagb) tau = 3.d0*(aj-tbagb)/(tn-tbagb)
         CALL star(kwp,mc,mc,tm,tn,tscls,lums,GB,zpars)
         lx = lmcgbf(mcx,GB)
         if(tau.lt.1.d0) lx = lums(2)*(lx/lums(2))**tau
         rx = rzhef(mc)
         rx = MIN(rhehgf(mc,lx,rx,lums(2)),rhegbf(lx))
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
         rc = rx
      elseif(kw.le.9)then
         if(wdflag.eq.0)then
            lx = 635.d0*mc*zpars(14)/((aco*0.1d0)**1.4d0)
         elseif(wdflag.ge.1)then
            lx = 300.d0*mc*zpars(14)/((aco*0.1d0)**1.18d0)
         endif
         rx = 0.0115d0*SQRT(MAX(1.48204d-06,
     &        (mch/mc)**(2.d0/3.d0) - (mc/mch)**(2.d0/3.d0)))
         rc = 5.d0*rx
      else
         rc = r
         menv = 1.0d-10
         renv = 1.0d-10
         k2 = 0.21d0
      endif
*
* Perturb the luminosity and radius due to small envelope mass.
*
      if(kw.ge.2.and.kw.le.9.and.kw.ne.7)then
         mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
         if(kw.ge.8) mew = ((mtc-mc)/mtc)*5.d0
         if(mew.lt.1.d0)then
            xx = lpertf(mt,mew)
            lum = lx*(lum/lx)**xx
            if(r.le.rx)then
               xx = 0.d0
            else
               xx = rpertf(mt,mew,r,rx)
            endif
            r = rx*(r/rx)**xx
         endif
         rc = MIN(rc,r)
      endif
*
* Calculate mass and radius of convective envelope, and envelope
* gyration radius.
*
      if(kw.lt.10)then
         CALL mrenv(kw,mass,mt,mc,lum,r,rc,aj,tm,lums(2),lums(3),
     &              lums(4),rzams,rtms,rg,menv,renv,k2)
      endif
*
      if(mass.gt.99.99d0)then
         mass = mass0
      endif
      if(mt.gt.99.99d0)then
         mt = mt0
      endif
*
      return
      end
***
