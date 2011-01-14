***
      SUBROUTINE star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
*
*
*       Stellar luminosity & evolution time. 
*       ------------------------------------
*
      implicit none
*
      integer kw
*
      real*8 mass,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
      real*8 tgb,tbagb,mch,mcmax,mc1,mc2,mcbagb,dx,am
      real*8 lambda,tau,mtc,mass0
      parameter(mch=1.44d0)
*
      real*8 lzamsf,lzahbf,lzhef
      real*8 tbgbf,thookf,tHef,themsf,mcgbf,mcagbf,mcheif,mcgbtf
      real*8 ltmsf,lbgbf,lHeIf,lHef,lbagbf,lmcgbf
      external lzamsf,lzahbf,lzhef
      external tbgbf,thookf,tHef,themsf,mcgbf,mcagbf,mcheif,mcgbtf
      external ltmsf,lbgbf,lHeIf,lHef,lbagbf,lmcgbf
*
*       Computes the characteristic luminosities at different stages (LUMS),
*       and various timescales (TSCLS).
*       Ref: P.P. Eggleton, M.J. Fitchett & C.A. Tout (1989) Ap.J. 347, 998.
*
*       Revised 27th March 1995 by C. A. Tout
*       and 24th October 1995 to include metallicity
*       and 13th December 1996 to include naked helium stars
*
*       Revised 5th April 1997 by J. R. Hurley
*       to include Z=0.001 as well as Z=0.02, convective overshooting,
*       MS hook and more elaborate CHeB. It now also sets the Giant
*       Branch parameters relevant to the mass of the star.
*
*       ------------------------------------------------------------
*       Times: 1; BGB              2; He ignition   3; He burning
*              4; Giant t(inf1)    5; Giant t(inf2) 6; Giant t(Mx)
*              7; FAGB t(inf1)     8; FAGB t(inf2)  9; FAGB  t(Mx)
*             10; SAGB t(inf1)    11; SAGB t(inf2) 12; SAGB  t(Mx)
*             13; TP              14; t(Mcmax)     
*
*       LUMS:  1; ZAMS             2; End MS        3; BGB
*              4; He ignition      5; He burning    6; L(Mx)
*              7; BAGB             8; TP
*
*       GB:    1; effective A(H)   2; A(H,He)       3; B
*              4; D                5; p             6; q
*              7; Mx               8; A(He)         9; Mc,BGB
*
*       ------------------------------------------------------------
*
*
      mass0 = mass
      if(mass0.gt.100.d0) mass = 100.d0
*
      if(kw.ge.7.and.kw.le.9) goto 90
      if(kw.ge.10) goto 95
*
* MS and BGB times
*
      tscls(1) = tbgbf(mass)
      tm = MAX(zpars(8),thookf(mass))*tscls(1)
*
* Zero- and terminal age main sequence luminosity
*
      lums(1) = lzamsf(mass)
      lums(2) = ltmsf(mass)
*
* Set the GB parameters
*
      GB(1) = MAX(-4.8d0,MIN(-5.7d0+0.8d0*mass,-4.1d0+0.14d0*mass))
      GB(1) = 10.d0**GB(1)
      GB(2) = 1.27d-05
      GB(8) = 8.0d-05
      GB(3) = MAX(3.0d+04,500.d0 + 1.75d+04*mass**0.6d0)
      if(mass.le.2.0)then
         GB(4) = zpars(6)
         GB(5) = 6.d0
         GB(6) = 3.d0
      elseif(mass.lt.2.5)then
         dx = zpars(6) - (0.975d0*zpars(6) - 0.18d0*2.5d0)
         GB(4) = zpars(6) - dx*(mass - 2.d0)/(0.5d0)
         GB(5) = 6.d0 - (mass - 2.d0)/(0.5d0)
         GB(6) = 3.d0 - (mass - 2.d0)/(0.5d0)
      else
         GB(4) = MAX(-1.d0,0.5d0*zpars(6) - 0.06d0*mass)
         GB(4) = MAX(GB(4),0.975d0*zpars(6) - 0.18d0*mass)
         GB(5) = 5.d0
         GB(6) = 2.d0
      endif
      GB(4) = 10.d0**GB(4)
      GB(7) = (GB(3)/GB(4))**(1.d0/(GB(5)-GB(6)))
*
* Change in slope of giant L-Mc relation.
      lums(6) = GB(4)*GB(7)**GB(5)
*
* HeI ignition luminosity
      lums(4) = lHeIf(mass,zpars(2)) 
      lums(7) = lbagbf(mass,zpars(2))
*
      if(mass.lt.0.1d0.and.kw.le.1)then
         tscls(2) = 1.1d0*tscls(1)
         tscls(3) = 0.1d0*tscls(1)
         lums(3) = lbgbf(mass) 
         goto 96
      endif
*
      if(mass.le.zpars(3))then
* Base of the giant branch luminosity
         lums(3) = lbgbf(mass) 
* Set GB timescales 
         tscls(4) = tscls(1) + (1.d0/((GB(5)-1.d0)*GB(1)*GB(4)))*
     &              ((GB(4)/lums(3))**((GB(5)-1.d0)/GB(5)))
         tscls(6) = tscls(4) - (tscls(4) - tscls(1))*((lums(3)/lums(6))
     &              **((GB(5)-1.d0)/GB(5)))
         tscls(5) = tscls(6) + (1.d0/((GB(6)-1.d0)*GB(1)*GB(3)))*
     &              ((GB(3)/lums(6))**((GB(6)-1.d0)/GB(6)))
* Set Helium ignition time
         if(lums(4).le.lums(6))then
            tscls(2) = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(1)*GB(4)))*
     &                      ((GB(4)/lums(4))**((GB(5)-1.d0)/GB(5)))
         else
            tscls(2) = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(1)*GB(3)))*
     &                      ((GB(3)/lums(4))**((GB(6)-1.d0)/GB(6)))
         endif
         tgb = tscls(2) - tscls(1)
         if(mass.le.zpars(2))then
            mc1 = mcgbf(lums(4),GB,lums(6))
            mc2 = mcagbf(mass)
            lums(5) = lzahbf(mass,mc1,zpars(2))
            tscls(3) = tHef(mass,mc1,zpars(2))
         else
            lums(5) = lHef(mass)*lums(4)
            tscls(3) = tHef(mass,1.d0,zpars(2))*tscls(1)
         endif
      else
* Note that for M>zpars(3) there is no GB as the star goes from
* HG -> CHeB -> AGB. So in effect tscls(1) refers to the time of
* Helium ignition and not the BGB.
         tscls(2) = tscls(1)
         tscls(3) = tHef(mass,1.d0,zpars(2))*tscls(1)
* This now represents the luminosity at the end of CHeB, ie. BAGB
         lums(5) = lums(7)
* We set lums(3) to be the luminosity at the end of the HG
         lums(3) = lums(4)
      endif
*
* Set the core mass at the BGB.
*
      if(mass.le.zpars(2))then
         GB(9) = mcgbf(lums(3),GB,lums(6))
      elseif(mass.le.zpars(3))then
         GB(9) = mcheif(mass,zpars(2),zpars(9))
      else
         GB(9) = mcheif(mass,zpars(2),zpars(10))
      endif
*
* FAGB time parameters
*
      tbagb = tscls(2) + tscls(3)
      tscls(7) = tbagb + (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &               ((GB(4)/lums(7))**((GB(5)-1.d0)/GB(5)))
      tscls(9) = tscls(7) - (tscls(7) - tbagb)*((lums(7)/lums(6))
     &                                    **((GB(5)-1.d0)/GB(5)))
      tscls(8) = tscls(9) + (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &               ((GB(3)/lums(6))**((GB(6)-1.d0)/GB(6)))
*
* Now to find Ltp and ttp using Mc,He,tp
*
      mcbagb = mcagbf(mass)
      mc1 = mcbagb
      if(mc1.ge.0.8d0.and.mc1.lt.2.25d0)then
* The star undergoes dredge-up at Ltp causing a decrease in Mc,He
         mc1 = 0.44d0*mc1 + 0.448d0
      endif
      lums(8) = lmcgbf(mc1,GB)
      if(mc1.le.GB(7))then
         tscls(13) = tscls(7) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                   (mc1**(1.d0-GB(5)))
      else
         tscls(13) = tscls(8) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                   (mc1**(1.d0-GB(6)))
      endif
*
* SAGB time parameters
*
      if(mc1.le.GB(7))then
         tscls(10) = tscls(13) + (1.d0/((GB(5)-1.d0)*GB(2)*GB(4)))*
     &               ((GB(4)/lums(8))**((GB(5)-1.d0)/GB(5)))
         tscls(12) = tscls(10) - (tscls(10) - tscls(13))*
     &               ((lums(8)/lums(6))**((GB(5)-1.d0)/GB(5)))
         tscls(11) = tscls(12) + (1.d0/((GB(6)-1.d0)*GB(2)*GB(3)))*
     &               ((GB(3)/lums(6))**((GB(6)-1.d0)/GB(6)))
      else
         tscls(10) = tscls(7)
         tscls(12) = tscls(9)
         tscls(11) = tscls(13) + (1.d0/((GB(6)-1.d0)*GB(2)*GB(3)))*
     &               ((GB(3)/lums(8))**((GB(6)-1.d0)/GB(6)))
      endif
*
* Get an idea of when Mc,C = Mc,C,max on the AGB
      tau = tscls(2) + tscls(3)
      mc2 = mcgbtf(tau,GB(8),GB,tscls(7),tscls(8),tscls(9))
      mcmax = MAX(MAX(mch,0.773d0*mcbagb - 0.35d0),1.05d0*mc2)
*
      if(mcmax.le.mc1)then
         if(mcmax.le.GB(7))then
            tscls(14) = tscls(7) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                      (mcmax**(1.d0-GB(5)))
         else
            tscls(14) = tscls(8) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                      (mcmax**(1.d0-GB(6)))
         endif
      else
* Star is on SAGB and we need to increase mcmax if any 3rd
* dredge-up has occurred.
         lambda = MIN(0.9d0,0.3d0+0.001d0*mass**5)
         mcmax = (mcmax - lambda*mc1)/(1.d0 - lambda)
         if(mcmax.le.GB(7))then
            tscls(14) = tscls(10) - (1.d0/((GB(5)-1.d0)*GB(2)*GB(4)))*
     &                      (mcmax**(1.d0-GB(5)))
         else
            tscls(14) = tscls(11) - (1.d0/((GB(6)-1.d0)*GB(2)*GB(3)))*
     &                      (mcmax**(1.d0-GB(6)))
         endif
      endif
      tscls(14) = MAX(tbagb,tscls(14))
      if(mass.ge.100.d0)then
         tn = tscls(2)
         goto 100
      endif
*
* Calculate the nuclear timescale - the time of exhausting
* nuclear fuel without further mass loss.
* This means we want to find when Mc = Mt which defines Tn and will
* be used in determining the timestep required. Note that after some 
* stars reach Mc = Mt there will be a Naked Helium Star lifetime
* which is also a nuclear burning period but is not included in Tn.
*
      if(ABS(mt-mcbagb).lt.1.0d-14.and.kw.lt.5)then
         tn = tbagb
      else
* Note that the only occurence of Mc being double-valued is for stars
* that have a dredge-up. If Mt = Mc where Mc could be the value taken
* from CHeB or from the AGB we need to check the current stellar type.
         if(mt.gt.mcbagb.or.(mt.ge.mc1.and.kw.gt.4))then
            if(kw.eq.6)then
               lambda = MIN(0.9d0,0.3d0+0.001d0*mass**5)
               mc1 = (mt - lambda*mc1)/(1.d0 - lambda)
            else
               mc1 = mt
            endif
            if(mc1.le.GB(7))then
               tn = tscls(10) - (1.d0/((GB(5)-1.d0)*GB(2)*GB(4)))*
     &                         (mc1**(1.d0-GB(5)))
            else
               tn = tscls(11) - (1.d0/((GB(6)-1.d0)*GB(2)*GB(3)))*
     &                         (mc1**(1.d0-GB(6)))
            endif
         else
            if(mass.gt.zpars(3))then
               mc1 = mcheif(mass,zpars(2),zpars(10))
               if(mt.le.mc1)then
                  tn = tscls(2)
               else
                  tn = tscls(2) + tscls(3)*((mt - mc1)/(mcbagb - mc1))
               endif
            elseif(mass.le.zpars(2))then
               mc1 = mcgbf(lums(3),GB,lums(6))
               mc2 = mcgbf(lums(4),GB,lums(6))
               if(mt.le.mc1)then
                  tn = tscls(1)
               elseif(mt.le.mc2)then
                  if(mt.le.GB(7))then
                     tn = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(1)*GB(4)))*
     &                               (mt**(1.d0-GB(5)))
                  else
                     tn = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(1)*GB(3)))*
     &                               (mt**(1.d0-GB(6)))
                  endif
               else
                  tn = tscls(2) + tscls(3)*((mt - mc2)/(mcbagb - mc2))
               endif
            else
               mc1 = mcheif(mass,zpars(2),zpars(9))
               mc2 = mcheif(mass,zpars(2),zpars(10))
               if(mt.le.mc1)then
                  tn = tscls(1)
               elseif(mt.le.mc2)then
                  tn = tscls(1) + tgb*((mt - mc1)/(mc2 - mc1))
               else
                  tn = tscls(2) + tscls(3)*((mt - mc2)/(mcbagb - mc2))
               endif
            endif
         endif
      endif
      tn = MIN(tn,tscls(14))
*
      goto 100
*
 90   continue
*
* Calculate Helium star Main Sequence lifetime.
*
      tm = themsf(mass)
      tscls(1) = tm
*
* Zero- and terminal age Helium star main sequence luminosity
*
      lums(1) = lzhef(mass)
      am = MAX(0.d0,0.85d0-0.08d0*mass)
      lums(2) = lums(1)*(1.d0+0.45d0+am)
*
* Set the Helium star GB parameters
*
      GB(8) = 8.0d-05
      GB(3) = 4.1d+04
      GB(4) = 5.5d+04/(1.d0+0.4d0*mass**4)
      GB(5) = 5.d0
      GB(6) = 3.d0
      GB(7) = (GB(3)/GB(4))**(1.d0/(GB(5)-GB(6)))
* Change in slope of giant L-Mc relation.
      lums(6) = GB(4)*GB(7)**GB(5)
*
*** Set Helium star GB timescales 
*
      mc1 = mcgbf(lums(2),GB,lums(6))
      tscls(4) = tm + (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                                          mc1**(1.d0-GB(5))
      tscls(6) = tscls(4) - (tscls(4) - tm)*((GB(7)/mc1)
     &                                     **(1.d0-GB(5)))
      tscls(5) = tscls(6) + (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                                       GB(7)**(1.d0-GB(6))
*
* Get an idea of when Mc = MIN(Mt,Mc,C,max) on the GB
      mtc = MIN(mt,1.45d0*mt-0.31d0)
      if(mtc.le.0.d0) mtc = mt
      mcmax = MIN(mtc,MAX(mch,0.773d0*mass-0.35d0))
      if(mcmax.le.GB(7))then
         tscls(14) = tscls(4) - (1.d0/((GB(5)-1.d0)*GB(8)*GB(4)))*
     &                   (mcmax**(1.d0-GB(5)))
      else
         tscls(14) = tscls(5) - (1.d0/((GB(6)-1.d0)*GB(8)*GB(3)))*
     &                   (mcmax**(1.d0-GB(6)))
      endif
      tscls(14) = MAX(tscls(14),tm)
      tn = tscls(14)
*
      goto 100
*
 95   continue
      tm = 1.0d+10
      tscls(1) = tm
 96   continue
      tn = 1.0d+10
*
 100  continue
      mass = mass0
*
      return
      end
***
