      SUBROUTINE TRDOT(I,DTM,M1)
*
*
*       Time-scale for expansion of radius.
*       -----------------------------------
*
      INCLUDE 'common6.h'
      REAL*8 TSCLS(20),LUMS(10),GB(10),TM,TN
      REAL*8 M0,M1,RM,LUM,AGE,MC,MC1,RCC,RM0,AGE0,M10
      REAL*8 menv,renv,k2
      REAL*8 pts1,pts2,eps,alpha2,tol
      PARAMETER(pts1=0.05d0,pts2=0.02d0)
      PARAMETER(eps=1.0d-06,alpha2=0.09d0,tol=1.0d-10)
*
*       Obtain stellar parameters at current epoch (body #I may be ghost).
      KW = KSTAR(I)
      M0 = BODY0(I)*ZMBAR
      IF(M1.LE.0.0) M1 = RADIUS(I)*SU
      M10 = M1
      MC = 0.D0
      AGE = TEV0(I)*TSTAR - EPOCH(I)
      CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
      CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &            RM,LUM,KW,MC,RCC,MENV,RENV,K2)
*
*       Quit if there is a change of type at the current TEV.
      if((kstar(i).le.6.and.kw.gt.6).or.
     &                            (kstar(i).le.9.and.kw.gt.9))then
         m1 = m10
         dtm = 0.d0
         goto 10
      endif
***
*
*       Base new time scale for changes in radius & mass on stellar type.
      if(kw.le.1)then
         dtm = pts1*tm
         dtr = tm - age
      elseif(kw.ge.10)then
         dtm = 1.0d+02
         dtr = dtm
      elseif(kw.eq.2)then
         dtm = pts1*(tscls(1) - tm)
         dtr = tscls(1) - age
      elseif(kw.eq.3)then
         if(age.lt.tscls(6))then
            dtm = pts2*(tscls(4) - age)
         else
            dtm = pts2*(tscls(5) - age)
         endif
         dtr = MIN(tscls(2),tn) - age
      elseif(kw.eq.4)then
         dtm = pts1*tscls(3)
         dtr = MIN(tn,tscls(2) + tscls(3)) - age
      elseif(kw.eq.5)then
         if(age.lt.tscls(9))then
            dtm = pts2*(tscls(7) - age)
         else
            dtm = pts2*(tscls(8) - age)
         endif
         dtr = MIN(tn,tscls(13)) - age
      elseif(kw.eq.6)then
         if(age.lt.tscls(12))then
            dtm = pts2*(tscls(10) - age)
         else
            dtm = pts2*(tscls(11) - age)
         endif
         dtm = MIN(dtm,0.005d0)
         dtr = tn - age
      elseif(kw.eq.7)then
         dtm = pts1*tm
         dtr = tm - age
      elseif(kw.eq.8.or.kw.eq.9)then
         if(age.lt.tscls(6))then
            dtm = pts2*(tscls(4) - age)
         else
            dtm = pts2*(tscls(5) - age)
         endif
         dtr = tn - age
      endif
*
* Record radius.
*
      rm0 = rm
      if(kw.ge.10) goto 30
      age0 = age
      kw0 = kw
      mc1 = mc
*
* Check for type change.
*
      it = 0
      if((dtr-dtm).le.tol)then
*
* Check final radius for too large a jump.
*
         age = MAX(age,age*(1.d0-eps) + dtr)
         CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RM,LUM,KW,MC1,RCC,MENV,RENV,K2)
         dr = rm - rm0
         if(ABS(dr).gt.0.1*rm0)then
            dtm = dtr - age0*eps
            dtdr = dtm/ABS(dr)
            dtm = alpha2*MAX(rm,rm0)*dtdr
            goto 20
         else
            dtm = dtr
            goto 30
         endif
      endif
*
* Limit to a 10% increase assuming no further mass loss
* and thus that the pertubation functions due to small envelope mass
* will not change the radius.
*
   20 age = age0 + dtm
      mc1 = mc
      CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &            RM,LUM,KW,MC1,RCC,MENV,RENV,K2)
      dr = rm - rm0
      it = it + 1
      if(it.eq.20.and.kw.eq.4) goto 30
      IF(IT.GT.30)THEN
         if(rank.eq.0) WRITE (6,22) IT, KSTAR(I), M0, DR, RM0
   22    FORMAT (' DANGER!    TRDOT: IT K* M0 DR RM0 ',2I4,1P,3E10.2)
         goto 30
      ENDIF
      if(ABS(dr).gt.0.1*rm0)then
         dtdr = dtm/ABS(dr)
         dtm = alpha2*MAX(rm0,rm)*dtdr
         if(it.ge.20) dtm = 0.5d0*dtm
         goto 20
      endif
*
 30   continue
*
*       Impose a lower limit and convert time interval to scaled units.
      DTM = MAX(DTM,1.0D-04)/TSTAR
*
*       Use dummy for routine CHINIT (new chain or chain collision).
   10 IF(IPHASE.EQ.8.OR.IPHASE.EQ.9)THEN
         M1 = RM
      ENDIF
*
      RETURN
      END
