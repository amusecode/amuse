***
      SUBROUTINE BSETID(kw,mass,mc,menv,r,rc,renv,lum,ospin,k2,q,
     &                  sep,ecc,oorb,delet1,dspin,eqspin,djt)
*
      implicit none
      integer kw
      real*8 mass,mc,menv,r,rc,renv,lum,ospin,k2,q
      real*8 sep,ecc,oorb,delet1,dspin,eqspin,djt
      real*8 renvk,rg2,k3
      PARAMETER(K3=0.21D0)
      real*8 ecc2,omecc2,sqome2,sqome3,raa2,raa6,twopi
      real*8 f1,f2,f3,f4,f5,tc,fc,tcqr,ttid
*
      ecc2 = ecc*ecc
      omecc2 = 1.d0 - ecc2
      sqome2 = SQRT(omecc2)
      sqome3 = sqome2**3
*
      raa2 = (r/sep)**2
      raa6 = raa2**3
      twopi = 2.d0*ACOS(-1.d0)
*
      f5 = 1.d0+ecc2*(3.d0+ecc2*0.375d0)
      f4 = 1.d0+ecc2*(1.5d0+ecc2*0.125d0)
      f3 = 1.d0+ecc2*(3.75d0+ecc2*(1.875d0+ecc2*7.8125d-02))
      f2 = 1.d0+ecc2*(7.5d0+ecc2*(5.625d0+ecc2*0.3125d0))
      f1 = 1.d0+ecc2*(15.5d0+ecc2*(31.875d0+ecc2*(11.5625d0
     &     +ecc2*0.390625d0)))
*
      if((kw.eq.1.and.mass.ge.1.25).or.kw.eq.4.or.kw.eq.7)then
         tc = 1.592d-09*(mass**2.84)
         fc = 1.9782d+04*SQRT((mass*r*r)/sep**5)
     &        *tc*(1.d0+q)**(5.d0/6.d0)
         tcqr = fc*q*raa6
         rg2 = k2
      elseif(kw.le.9)then
         renvk = MIN(renv,r-rc)
         renvk = MAX(renvk,1.0d-10)
         tc = (menv*renvk*(r-0.5d0*renvk)/
     &        (3.d0*lum))**(1.d0/3.d0)
         tc = 0.4311D0*tc
         ttid = twopi/(1.0d-10 + ABS(oorb - ospin))
         fc = MIN(1.d0,(ttid/(2.d0*tc)**2))
         tcqr = 2.d0*fc*q*raa6*menv/(21.d0*tc*mass)
         rg2 = (k2*(mass-mc))/mass
      else
         fc = 7.33d-09*(lum/mass)**(5.d0/7.d0)
         tcqr = fc*q*q*raa2*raa2/(1.d0+q)
         rg2 = k3
      endif
*
      delet1 = 27.d0*tcqr*(1.d0+q)*raa2*(ecc/sqome2**13)*
     &         (f3 - (11.d0/18.d0)*sqome3*f4*ospin/oorb)
      dspin = (3.d0*q*tcqr/(rg2*omecc2**6))*
     &        (f2*oorb - sqome3*f5*ospin)
      eqspin = oorb*f2/(sqome3*f5)
      djt = (k2*(mass-mc)*r*r + k3*rc*rc*mc)*dspin
*
      RETURN
      END
***
