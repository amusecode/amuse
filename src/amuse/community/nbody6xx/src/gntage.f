      SUBROUTINE gntage(mc,mt,kw,zpars,m0,aj)
*
* A routine to determine the age of a giant from its core mass and type.
* ----------------------------------------------------------------------
*
*       Author : C. A. Tout
*       Date   : 24th September 1996
*       Revised: 21st February 1997 to include core-helium-burning stars
*
*       Rewritten: 2nd January 1998 by J. R. Hurley to be compatible with
*                  new evolution routines and to include new stellar types.
*
      implicit none
*
      integer kw
      integer j,jmax
      parameter(jmax=30)
*
      real*8 mc,mt,m0,aj,tm,tn
      real*8 tscls(20),lums(10),GB(10),zpars(20)
      real*8 mmin,mmax,mmid,dm,f,fmid,dell,derl,lum
      real*8 macc,lacc,tiny
      parameter(macc=0.00001d0,lacc=0.0001d0,tiny=1.0d-14)
      real*8 mcx,mcy
*
      real*8 mcheif,mcagbf,mheif,mbagbf,mcgbf,lmcgbf,lbgbf,lbgbdf
      external mcheif,mcagbf,mheif,mbagbf,mcgbf,lmcgbf,lbgbf,lbgbdf
*
* This should only be entered with KW = 3, 4, 5, 6 or 9
*
* First we check that we don't have a CheB star 
* with too small a core mass.
      if(kw.eq.4)then
* Set the minimum CHeB core mass using M = Mflash
         mcy = mcheif(zpars(2),zpars(2),zpars(10))
         if(mc.le.mcy) kw = 3
      endif
*
* Next we check that we don't have a GB star for M => Mfgb
      if(kw.eq.3)then
* Set the maximum GB core mass using M = Mfgb
         mcy = mcheif(zpars(3),zpars(2),zpars(9))
         if(mc.ge.mcy)then
            kw = 4
            aj = 0.d0
         endif
      endif
*
      if(kw.eq.6)then
*
* We try to start the star from the start of the SAGB by
* setting Mc = Mc,TP.
*
         mcy = 0.44d0*2.25d0 + 0.448d0
         if(mc.gt.mcy)then
* A type 6 with this sized core mass cannot exist as it should
* already have become a NS or BH as a type 5. 
* We set it up so that it will.
            mcx = (mc + 0.35d0)/0.773d0
         elseif(mc.ge.0.8d0)then
            mcx = (mc - 0.448d0)/0.44d0
         else
            mcx = mc
         endif
         m0 = mbagbf(mcx)
         if(m0.lt.tiny)then
* Carbon core mass is less then the minimum for the start of SAGB.
* This must be the case of a low-mass C/O or O/Ne WD with only a
* very small envelope added or possibly the merger of a helium star
* with a main sequence star. We will set m0 = mt and then reset the
* core mass to allow for some helium to be added to the C/O core. 
            kw = 14
         else
            CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
            aj = tscls(13)
         endif
      endif
*
      if(kw.eq.5)then
*
* We fit a Helium core mass at the base of the AGB.
*
         m0 = mbagbf(mc)
         if(m0.lt.tiny)then
* Helium core mass is less then the BAGB minimum.
            kw = 14
         else
            CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
            aj = tscls(2) + tscls(3)
         endif
      endif
*
*
      if(kw.eq.4)then
*
* The supplied age is actually the fractional age, fage, of CHeB lifetime
* that has been completed, ie. 0 <= aj <= 1.
*
* Get the minimum, fage=1, and maximum, fage=0, allowable masses
         mcy = mcagbf(zpars(2))
         if(mc.ge.mcy)then
            mmin = mbagbf(mc)
         else
            mmin = zpars(2)
         endif
         mmax = mheif(mc,zpars(2),zpars(10))
         if(aj.lt.tiny)then
            m0 = mmax
            goto 20
         elseif(aj.ge.1.d0)then
            m0 = mmin
            goto 20
         endif
* Use the bisection method to find m0
         fmid = (1.d0-aj)*mcheif(mmax,zpars(2),zpars(10)) +
     &          aj*mcagbf(mmax) - mc
         f = (1.d0-aj)*mcheif(mmin,zpars(2),zpars(10)) +
     &       aj*mcagbf(mmin) - mc
         if(f*fmid.ge.0.d0)then
* This will probably occur if mc is just greater than the minimum
* allowed mass for a CHeB star and fage > 0.
            kw = 3
            goto 90
         endif
         m0 = mmin
         dm = mmax - mmin
         do 10 j = 1,jmax
            dm = 0.5d0*dm
            mmid = m0 + dm
            fmid = (1.d0-aj)*mcheif(mmid,zpars(2),zpars(10)) +
     &             aj*mcagbf(mmid) - mc
            if(fmid.lt.0.d0) m0 = mmid
            if(abs(dm).lt.macc.or.ABS(fmid).lt.tiny) goto 20
 10      continue
 20      continue
*
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         aj = tscls(2) + aj*tscls(3)
*
      endif
*
 90   continue
*
      if(kw.eq.3)then
*
* First we double check that we don't have a GB star for M => Mfgb
         mcy = mcheif(zpars(3),zpars(2),zpars(9))
* Next we find an m0 so as to place the star at the BGB
         mcx = mcheif(zpars(2),zpars(2),zpars(9))
         if(mc.gt.mcx)then
            m0 = mheif(mc,zpars(2),zpars(9))
         else
* Use Newton-Raphson to find m0 from Lbgb
            m0 = zpars(2)
            CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
            lum = lmcgbf(mc,GB)
            j = 0
 30         continue
            dell = lbgbf(m0) - lum
            if(ABS(dell/lum).le.lacc) goto 40
            derl = lbgbdf(m0)
            m0 = m0 - dell/derl
            j = j + 1
            goto 30
 40         continue
         endif
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         aj = tscls(1) + 1.0d-06*(tscls(2) - tscls(1))
*
      endif
*
      if(kw.eq.8.or.kw.eq.9)then
*
* We make a post-MS naked helium star.
* To make things easier we put the star at the TMS point
* so it actually begins as type 8.
*
         kw = 8
         mmin = mc
         CALL star(kw,mmin,mc,tm,tn,tscls,lums,GB,zpars)
         mcx = mcgbf(lums(2),GB,lums(6))
         f = mcx - mc
         mmax = mt
         do 50 j = 1,jmax
            CALL star(kw,mmax,mc,tm,tn,tscls,lums,GB,zpars)
            mcy = mcgbf(lums(2),GB,lums(6))
            if(mcy.gt.mc) goto 60
            mmax = 2.d0*mmax
 50      continue
 60      continue
         fmid = mcy - mc
* Use the bisection method to find m0
         m0 = mmin
         dm = mmax - mmin
         do 70 j = 1,jmax
            dm = 0.5d0*dm
            mmid = m0 + dm
            CALL star(kw,mmid,mc,tm,tn,tscls,lums,GB,zpars)
            mcy = mcgbf(lums(2),GB,lums(6))
            fmid = mcy - mc
            if(fmid.lt.0.d0) m0 = mmid
            if(abs(dm).lt.macc.or.ABS(fmid).lt.tiny) goto 80
 70      continue
 80      continue
*
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         aj = tm + 1.0d-10*tm
*
      endif
*
      if(kw.eq.14)then
*
         kw = 4
         m0 = mt
         mcy = mcagbf(m0)
         aj = mc/mcy
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         if(m0.le.zpars(2))then
            mcx = mcgbf(lums(4),GB,lums(6))
         else
            mcx = mcheif(m0,zpars(2),zpars(10))
         end if
         mc = mcx + (mcy - mcx)*aj
         aj = tscls(2) + aj*tscls(3)
      endif
*
      RETURN
      END
