      SUBROUTINE TRFLOW(IPAIR,DTR)
*
*
*       Time until Roche overflow.
*       --------------------------
*
*
      INCLUDE 'common6.h'
      REAL*8  TSCLS(20),LUMS(10),GB(10),TM,TN
      REAL*8  M0,M1,AGE,LUM,MC,RCC,Q,RL1,RL2
      REAL*8  MENV,RENV,K2
      REAL*8  AJ,RHIGH,THIGH,T,TT,RR,RX,TLOW,RLOW,DEL,DER,EPS,TOL
      PARAMETER (EPS = 1.0D-06, TOL = 1.0D-04)
      REAL*8 RL,RTMSF,RGBF,RGBDF,RAGBF,RAGBDF,RZHEF
      EXTERNAL RL,RTMSF,RGBF,RGBDF,RAGBF,RAGBDF,RZHEF
      LOGICAL ITEV
*
*
*       Set index of first component and semi-major axis.
      ITEV = .TRUE.
      IF(IPAIR.LT.0)THEN
         ITEV = .FALSE.
         IPAIR = -IPAIR
      ENDIF
      I1 = 2*IPAIR - 1
      I = N + IPAIR
      SEMI = -0.5*BODY(I)/H(IPAIR)
*
*       Determine indices for primary & secondary star (donor & accretor).
      J1 = 2*IPAIR - 1
      J2 = I1 + 1
      Q = BODY(J1)/BODY(J2)
      RL1 = RL(Q)*SEMI
*       Evaluate Roche radius for the second star.
      Q = 1.0/Q
      RL2 = RL(Q)*SEMI
*
*       Compare scaled Roche radii when choosing the primary.
      IF (RADIUS(J1)/RL1.LT.RADIUS(J2)/RL2) THEN
          J1 = I1 + 1
          J2 = I1
          RL1 = RL2
      END IF
*     IF(TEV(J1).LE.0.0)THEN
*        WRITE(6,99)J1,NAME(J1),NAME(I),KSTAR(J1),TEV(J1)+TOFF,TTOT
*99      FORMAT(' TRFLOW WARNING! ',3I6,I4,2F10.3)
*     ENDIF
*
*       Exit with large interval if semi-major axis is negative.
      IF(SEMI.LE.0.D0)THEN
          DTR = 1.0D+10
          GOTO 210
*       Assign DTR for active Roche during coasting.
      ELSEIF(KSTAR(I).GT.0.AND.MOD(KSTAR(I),2).NE.0)THEN
          IF(TEV(I).GT.9.9D+09) TEV(I) = TEV0(J1)
          DTR = TEV(I) - TIME
          GOTO 210
      ENDIF
*
*       Convert Roche radius, radius and initial & current mass to SU.
      RLS = RL1*SU
      M0 = BODY0(J1)*ZMBAR
      M1 = BODY(J1)*ZMBAR
      MC = 0.D0
      KW = KSTAR(J1)
      CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
      AGE = TEV0(J1)*TSTAR - EPOCH(J1)
      CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &            RSJ,LUM,KW,MC,RCC,MENV,RENV,K2)
*
* If the star already fills its Roche lobe exit with DTR = 0.
*
      IF(RSJ.GE.RLS)THEN
         DTR = 0.D0
         RADIUS(J1) = RSJ/SU
         GOTO 210
      ENDIF
*
*       Exit with large interval if primary has terminated its evolution.
      IF(KSTAR(J1).GE.10)THEN
          DTR = 1.0D+10
          GOTO 210
      ENDIF
*
* If the star is on the MS then see if it will fill its Roche Lobe
* during this stage.
*
      T = 2.0D+10
      IF(KW.LE.1.OR.KW.EQ.7)THEN
         AJ = 0.99D0*TM
         CALL hrdiag(M0,AJ,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RHIGH,LUM,KW,MC,RCC,MENV,RENV,K2)
         IF(RHIGH.GT.RLS.AND.AJ.GT.AGE)THEN
* Find the overflow time using binary chopping.
            TLOW = AGE
            RLOW = RSJ
            THIGH = AJ
            IT = 0
   10       T = 0.5D0*(THIGH + TLOW)
            CALL hrdiag(M0,T,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                  RR,LUM,KW,MC,RCC,MENV,RENV,K2)
            IF(RR.LT.RLS)THEN
               TLOW = T
               RLOW = RR
            ELSE
* Stop when 0.1% accuracy is achieved.
               IT = IT + 1
               IF(IT.EQ.25)THEN
                  IF(ABS(RR-RLS).GT.0.1)THEN
                     if(rank.eq.0)WRITE(38,*)' TRFLOW KW1: ',rr,rls,t,tm
                  ENDIF
               ENDIF
               IF(ABS(RR - RLS).LE.0.001*RR.OR.IT.GT.25) GOTO 200
               THIGH = T
               RHIGH = RR
            ENDIF
            GOTO 10
         ENDIF
      ENDIF
*
* See if the star will fill its Roche Lobe on the HG if not already
* past this stage.
*
      IF(KW.LE.2)THEN
         AJ = MIN(TN,TSCLS(1))*(1.D0-EPS)
         CALL hrdiag(M0,AJ,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RHIGH,LUM,KW,MC,RCC,MENV,RENV,K2)
         IF(RHIGH.GT.RLS.AND.AJ.GT.AGE)THEN
* Solve for the overflow time.
            RR = RTMSF(M0)
            T = LOG(RLS/RR)/LOG(RHIGH/RR)
            T = TM + T*(TSCLS(1) - TM)
            GOTO 200
         ENDIF
         IF(TN.LE.TSCLS(1)) GOTO 200
      ENDIF
*
* If the star is pre-helium ignition see if it will fill its Roche
* lobe before helium ignition.
*
      IF(KW.LE.3)THEN
         AJ = MIN(TN,TSCLS(2))*(1.D0-EPS)
         CALL hrdiag(M0,AJ,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RHIGH,LUM,KW,MC,RCC,MENV,RENV,K2)
         IF(RHIGH.GE.RLS)THEN
* Solve for the luminosity corresponding to the Roche radius and
* then for the overflow time.
            IT = 0
            LUM = LUMS(3)
 30         IT = IT + 1
            DEL = RGBF(M1,LUM) - RLS
            IF(IT.EQ.25.AND.ABS(DEL).GT.0.1)THEN
               if(rank.eq.0)then
               WRITE(38,*)' TRFLOW KW3: ',rgbf(m1,lum),rls,lum,lums(3)
               end if
            ENDIF
            IF(ABS(DEL/RLS).LE.TOL.OR.IT.EQ.25) GOTO 40
            DER = RGBDF(M1,LUM)
            LUM = LUM - DEL/DER
            GOTO 30
 40         CONTINUE
            IF(LUM.LE.LUMS(6))THEN
               T = TSCLS(4) - (1.D0/((GB(5)-1.D0)*GB(1)*GB(4)))*
     &                      ((GB(4)/LUM)**((GB(5)-1.D0)/GB(5)))
            ELSE
               T = TSCLS(5) - (1.D0/((GB(6)-1.D0)*GB(1)*GB(3)))*
     &                      ((GB(3)/LUM)**((GB(6)-1.D0)/GB(6)))
            ENDIF
            GOTO 200
         ELSE
* If a low mass star has not yet filled its Roche Lobe then it will
* not do it until after the Helium Flash. Thus we don't let it go
* any further until it has actually become a type 4.
            IF(M0.LT.ZPARS(2).OR.TN.LE.TSCLS(2)) GOTO 200
         ENDIF
      ENDIF
*
* Check for overflow during the CHeB stage.
*
      IF(KW.EQ.4)THEN
*
         IF(TN.LT.(TSCLS(2)+TSCLS(3)))THEN
            T = MAX(TSCLS(2),AGE)
            AJ = T + 0.5D0*(TN - T)
            CALL hrdiag(M0,AJ,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                  RHIGH,LUM,KW,MC,RCC,MENV,RENV,K2)
            IF(RHIGH.LT.RLS)THEN
* If the evolution is due to end during CHeB then the pertubation
* functions for small envelope mass will take effect at some point
* causing the radius to decrease. We assume that this point is after
* AJ and quit with T as an underestimate of the overflow time.
               T = AJ
               GOTO 200
            ENDIF
         ELSE
            AJ = (TSCLS(2)+TSCLS(3))*(1.D0-EPS)
            CALL hrdiag(M0,AJ,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                  RHIGH,LUM,KW,MC,RCC,MENV,RENV,K2)
            IF(RHIGH.LT.RLS) GOTO 50
         ENDIF
* Find the overflow time using binary chopping.
         TLOW = AGE
         RLOW = RSJ
         THIGH = AJ
         IT = 0
 60      T = 0.5D0*(THIGH + TLOW)
         CALL hrdiag(M0,T,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RR,LUM,KW,MC,RCC,MENV,RENV,K2)
         IF(RR.LT.RLS)THEN
            TLOW = T
            RLOW = RR
         ELSE
* Stop when 0.1% accuracy is achieved.
            IT = IT + 1
            IF(IT.EQ.25)THEN
               IF(rank.eq.0.and.ABS(RR-RLS).GT.0.1)THEN
                  WRITE(38,*)' TRFLOW KW4: ',rr,rls,tscls(2),t,tscls(3)
               ENDIF
            ENDIF
            IF(ABS(RR - RLS).LE.0.001*RR.OR.IT.GT.25) GOTO 200
            THIGH = T
            RHIGH = RR
         ENDIF
         GOTO 60
 50      CONTINUE
      ENDIF
*
* The star now has only until the end of the AGB phase to fill its lobe.
*
      IF(KW.LE.6)THEN
* Solve for the luminosity corresponding to the Roche radius and
* then for the overflow time.
         IT = 0
         IF(KW.EQ.6)THEN
            LUM = LUMS(8)
         ELSE
            LUM = LUMS(7)
         ENDIF
         RR = RAGBF(M1,LUM,ZPARS(2))
         IF(RR.GT.RLS)THEN
* In this case the solution should already have been found. If it hasn't
* then most likely the envelope is small and the pertubation functions
* are taking effect.
            T = TN
            GOTO 200
         ENDIF
 70      IT = IT + 1
         RR = RAGBF(M1,LUM,ZPARS(2))
         DEL = RR - RLS
         IF(IT.EQ.25.AND.ABS(DEL).GT.0.1)THEN
            if(rank.eq.0) WRITE(38,*)' TRFLOW KW6: ',rr,rls,lum,lums(7)
         ENDIF
         IF(ABS(DEL/RLS).LE.TOL.OR.IT.EQ.25) GOTO 80
         DER = RAGBDF(M1,LUM,ZPARS(2))
         LUM = LUM - DEL/DER
         GOTO 70
 80      CONTINUE
         IF(LUM.LE.LUMS(8))THEN
            IF(LUM.LE.LUMS(6))THEN
               T = TSCLS(7) - (1.D0/((GB(5)-1.D0)*GB(8)*GB(4)))*
     &                      ((GB(4)/LUM)**((GB(5)-1.D0)/GB(5)))
            ELSE
               T = TSCLS(8) - (1.D0/((GB(6)-1.D0)*GB(8)*GB(3)))*
     &                      ((GB(3)/LUM)**((GB(6)-1.D0)/GB(6)))
            ENDIF
         ELSE
            IF(LUM.LE.LUMS(6))THEN
               T = TSCLS(10) - (1.D0/((GB(5)-1.D0)*GB(2)*GB(4)))*
     &                       ((GB(4)/LUM)**((GB(5)-1.D0)/GB(5)))
            ELSE
               T = TSCLS(11) - (1.D0/((GB(6)-1.D0)*GB(2)*GB(3)))*
     &                       ((GB(3)/LUM)**((GB(6)-1.D0)/GB(6)))
            ENDIF
         ENDIF
         GOTO 200
      ENDIF
*
      IF(KW.LE.9)THEN
*
* The star has until the end of the Helium GB to fill its lobe.
*
         LUM = (RLS/0.08D0)**(4.0/3.0)
         CM = 2.0D-03*M1**2.5/(2.D0 + M1**5)
         IF(CM*LUM.GE.100.D0)THEN
            DTR = 1.0D+10
            GOTO 210
         ENDIF
         RX = RZHEF(M1)
         RR = RX*(LUM/LUMS(2))**0.2 + 
     &            0.02D0*(EXP(CM*LUM) - EXP(CM*LUMS(2)))
         IF(RR.LT.RLS)THEN
* Find the overflow luminosity using binary chopping.
            TLOW = LUM
            RLOW = RR
 89         LUM = 2.0*LUM
            RR = RX*(LUM/LUMS(2))**0.2 + 
     &               0.02D0*(EXP(CM*LUM) - EXP(CM*LUMS(2)))
            IF(RR.LE.RLS) GOTO 89
            THIGH = LUM
            RHIGH = RR
            IT = 0
 90         LUM = 0.5D0*(THIGH + TLOW)
            RR = RX*(LUM/LUMS(2))**0.2 + 
     &               0.02D0*(EXP(CM*LUM) - EXP(CM*LUMS(2)))
            IF(RR.LT.RLS)THEN
               TLOW = LUM
               RLOW = RR
            ELSE
* Stop when 0.1% accuracy is achieved.
               IT = IT + 1
               IF(IT.EQ.25)THEN
                  IF(ABS(RR-RLS).GT.0.1)THEN
                     if(rank.eq.0)WRITE(38,*)' TRFLOW KW9: ',rr,rls,t,tm
                  ENDIF
               ENDIF
               IF(ABS(RR - RLS).LE.0.001*RR.OR.IT.GT.25) GOTO 95
               THIGH = LUM
               RHIGH = RR
            ENDIF
            GOTO 90
 95         CONTINUE
         ENDIF
         IF(LUM.LE.LUMS(6))THEN
            T = TSCLS(4) - (1.D0/((GB(5)-1.D0)*GB(8)*GB(4)))*
     &                   ((GB(4)/LUM)**((GB(5)-1.D0)/GB(5)))
         ELSE
            T = TSCLS(5) - (1.D0/((GB(6)-1.D0)*GB(8)*GB(3)))*
     &                   ((GB(3)/LUM)**((GB(6)-1.D0)/GB(6)))
         ENDIF
      ENDIF
*
 200  CONTINUE
      IF(T.GE.TN)THEN
         DTR = 1.0D+10
      ELSE
         TT = (T+EPOCH(J1))/TSTAR
         DTR = TT - TIME
         IF(DTR.LT.0.D0.AND.ITEV)THEN
            TEV(J1) = MIN(TEV(J1),TIME)
            TEV(J2) = TEV(J1)
         ENDIF
         DTR = MAX(DTR,0.D0)
         IF(TEV(J1).LT.TEV0(J1))THEN
            if(rank.eq.0) WRITE (6,101)J1,KSTAR(J1),TEV(J1)+TOFF,TTOT
 101        FORMAT(' TRFLOW WARNING! TEV<TEV0 ',I6,I4,2F10.3)
         ENDIF
      ENDIF
*
 210  RETURN
      END
***
