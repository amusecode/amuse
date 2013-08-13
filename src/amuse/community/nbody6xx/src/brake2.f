      SUBROUTINE BRAKE2(IPAIR,ITERM)
*
*
*       Gravitational radiation of hierarchical binary.
*       -----------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/SLOW0/  RANGE,ISLOW(10)
      REAL*8 M0,AJ,R1,LUM1,MC,RCC,MENV,RENV,K2
      REAL*8 TSCLS(20),LUMS(10),GB(10),TM,TN
      REAL*8  M1,M2,RL,Q
      EXTERNAL  RL,RTMSF
*
*
      I = N + IPAIR
      I1 = 2*IPAIR - 1
      IM = 0
      DO 1 K = 1,NMERGE
         IF (NAMEM(K).EQ.NAME(I)) IM = K
    1 CONTINUE
*       Find location of the secondary (i.e. ghost).
      CALL FINDJ(I,I2,IM)
      M1 = CM(1,IM)*ZMBAR
      M2 = CM(2,IM)*ZMBAR
*
*       Quit if active ROCHE or mass-losing star.
      ITERM = 1
      IF(KSTARM(IM).EQ.11.OR.KSTARM(IM).EQ.13)THEN
          GO TO 30
      END IF
      IF (M1.GE.10.0.OR.(KSTAR(I1).GT.1.AND.KSTAR(I1).LT.10))THEN
          GO TO 30
      END IF
      IF (M2.GE.10.0.OR.(KSTAR(I2).GT.1.AND.KSTAR(I2).LT.10))THEN
          GO TO 30
      END IF
*
*       Skip braking if semi-major axis exceeds 10 solar radii.
      ITERM = 0
      ZMB = CM(1,IM) + CM(2,IM)
      SEMI = -0.5*ZMB/HM(IM)
      IF (SEMI*SU.GT.10.0)THEN
          TGR = 1.0D+06
          TBR = 1.0D+06
          SEMI1 = SEMI
          GO TO 20
      END IF
*
*       Obtain current separation, square velocity and t'' = 2*U*U'.
      RB = 0.0
      V20 = 0.0
      TD2 = 0.0
      DO 5 K = 1,4
          RB = RB + UM(K,IM)**2
          V20 = V20 + UMDOT(K,IM)**2
          TD2 = TD2 + 2.0*UM(K,IM)*UMDOT(K,IM)
    5 CONTINUE
*
*       Evaluate inner eccentricity and outer period.
      ECC2 = (1.0 - RB/SEMI)**2 + TD2**2/(ZMB*SEMI)
      ECC = SQRT(ECC2)
      SEMI0 = -0.5*BODY(I)/H(IPAIR)
*
*       Identify the most magnetically active star (#J2).
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
          IF ((KSTAR(I1).EQ.1.AND.M1.GT.1.25).OR.
     &        (KSTAR(I1).EQ.0.AND.M1.LT.0.35)) THEN
              J2 = I2
          ELSE
              J2 = I1
          END IF
      ELSE
          IF ((KSTAR(I2).EQ.1.AND.M2.GT.1.25).OR.
     &        (KSTAR(I2).EQ.0.AND.M2.LT.0.35)) THEN
              J2 = I1
          ELSE
              J2 = I2
          END IF
      END IF
*
*       Set mass and radius of star #J2 in solar units.
      IF (J2.EQ.I1) M2 = M1
      R2 = RADIUS(J2)*SU
*
*       Define binary period, solar radius & GMS in cgs units.
      TB = 3.147D+07*YRS*SEMI*SQRT(SEMI/ZMB)
      RSUN = 6.96D+10
      GMS = GM*ZMBAR*ZMB
*
*       Form time derivative of angular momentum (Regos & Tout 1995).
      ZJDOT = 3.8D-30*1.989D+33*M2*RSUN**4*R2**3*(TWOPI/TB)**3
*
*       Determine new semi-major axis from angular momentum relation.
      ZMU = CM(1,IM)*CM(2,IM)/ZMB
      ACM = 3.08D+18*RBAR*SEMI
      ADOT = 2.0*SQRT(ACM/GMS)*ZJDOT/(1.989D+33*ZMBAR*ZMU)
*
*       Define old primary and evaluate time scale for spin-down (yrs).
      TBR = ACM/(3.147D+07*ADOT)
*       Evaluate timescale for 10% change in semi-major axis (N-body units). 
      TBR = 0.1D0*TBR/(1.0D+06*TSTAR)
*
*       Include gravitational radiation (64*GM**3/5*c**5; Ap.J. 418, 147).
*     FAC = 64.0*((6.67D-08*1.989D+33)/3.0D+10)**3/(5.0*9.0D+20)
      AGDOT = 1.23D+27*CM(1,IM)*CM(2,IM)*ZMB*(ZMBAR/ACM)**3
      TGR = ACM/(3.147D+07*AGDOT)
      TGR = 0.1D0*TGR/(1.0D+06*TSTAR)
*
*       Suppress magnetic braking for massive MS/low-mass or evolved star.
      IF (((M2.GT.1.25.OR.M2.LT.0.35).AND.KSTAR(J2).LE.2).OR.
     &      KSTAR(J2).GE.10) THEN
          ADOT = 0.0
          TBR = 1.0D+06
      END IF
*
*       Specify time interval.
      DT = TEV(I) - TEV0(I)
      TK = TWOPI*ABS(SEMI0)*SQRT(ABS(SEMI0)/BODY(I))
*
*       Combine effects but include possible cutoff of magnetic braking.
      ADOT = ADOT + AGDOT
*
*       Convert from cgs to scaled units and update semi-major axis.
      ADOT = ADOT/(1.0D+05*VSTAR)
      SEMI1 = SEMI - ADOT*DT
*
*       Include safety test on new semi-major axis.
      IF (ABS(SEMI1).LT.RADIUS(J2).OR.SEMI1.LT.0.0) THEN
          SEMI1 = RADIUS(J2)
      END IF
*
*       Update binding energy.
      HI = HM(IM)
      HM(IM) = -0.5*ZMB/SEMI1
*
*       Correct EMERGE & ECOLL (consistency; no net effect).
      DECORR = ZMU*(HI - HM(IM))
      EMERGE = EMERGE - DECORR
      ECOLL = ECOLL + DECORR
      EGRAV = EGRAV + DECORR
*
*       Specify KS coordinate & velocity scaling factors at general point.
      C2 = SQRT(SEMI1/SEMI)
      V2 = 0.5*(ZMB + HM(IM)*RB*(SEMI1/SEMI))
      C1 = SQRT(V2/V20)
*
*       Re-scale KS variables to new energy with constant eccentricity.
      RB = 0.0D0
      DO 10 K = 1,4
          UM(K,IM) = C2*UM(K,IM)
          UMDOT(K,IM) = C1*UMDOT(K,IM)
          RB = RB + UM(K,IM)**2
   10 CONTINUE
*
*       Rectify the hierarchical KS variables.
      CALL HIRECT(IM)
*
*       Determine indices for primary & secondary star (donor & accretor).
   20 J1 = I1
      J2 = I2
      Q = CM(1,IM)/CM(2,IM)
      RL1 = RL(Q)*SEMI1
*       Evaluate Roche radius for the second star.
      Q = 1.0/Q
      RL2 = RL(Q)*SEMI1
*
*       Update radius for MS stars to current time. 
      TM0 = 100.D0
      RTMS = 0.D0
      IF(KSTAR(I1).LE.1)THEN
         M0 = M1
         MC = 0.D0
         AJ = TIME*TSTAR - EPOCH(I1)
         KW = KSTAR(I1)
         CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
         CALL hrdiag(M0,AJ,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               R1,LUM1,KW,MC,RCC,MENV,RENV,K2)
         IF(KW.GT.1)THEN
            TM0 = 0.D0
            TEV(I1) = MIN(TEV(I1),TIME)
         ELSE
            TM = (TM + EPOCH(I1))/TSTAR - TIME
            TM0 = MIN(TM0,TM)
            RADIUS(I1) = R1/SU
            TEV0(I1) = TIME
            RTMS = RTMSF(M1)/SU
            DTMS = TM
         ENDIF
      ENDIF
      IF(KSTAR(I2).LE.1)THEN
         M0 = M2
         MC = 0.D0
         AJ = TIME*TSTAR - EPOCH(I2)
         KW = KSTAR(I2)
         CALL star(KW,M0,M2,TM,TN,TSCLS,LUMS,GB,ZPARS)
         CALL hrdiag(M0,AJ,M2,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               R1,LUM1,KW,MC,RCC,MENV,RENV,K2)
         IF(KW.GT.1)THEN
            TM0 = 0.D0
            TEV(I2) = MIN(TEV(I2),TIME)
         ELSE
            TM = (TM + EPOCH(I2))/TSTAR - TIME
            TM0 = MIN(TM0,TM)
            RADIUS(I2) = R1/SU
            TEV0(I2) = TIME
            RTMS2 = RTMSF(M2)/SU
            IF(RTMS2.GT.RTMS)THEN
               RTMS = RTMS2
               DTMS = TM
            ENDIF
         ENDIF
      ENDIF
      IF(KSTAR(I1).GE.10.AND.TM0.GT.0.D0) TEV0(I1) = TIME
      IF(KSTAR(I2).GE.10.AND.TM0.GT.0.D0) TEV0(I2) = TIME
*
*       Compare scaled Roche radii when choosing the primary.
      IF (RADIUS(J1)/RL1.LT.RADIUS(J2)/RL2) THEN
          J1 = I2
          J2 = I1
          RL1 = RL2
      END IF
*
*       Check for Roche mass transfer at every stage (AM CVn formation).
      IF (RADIUS(J1).GT.RL1.AND.KZ(34).GT.0) THEN
          WRITE (6,25)  NAME(I1), NAMEG(IM), KSTAR(I1), KSTAR(I2),
     &                  SEMI1*SU, RL1*SU, RADIUS(J1)*SU
   25     FORMAT (' BRAKE2 TERM    NM K* A RL R*  ',2I6,2I4,3F7.3)
          ITERM = 1
*
*       Check for termination of MS stage. 
      ELSEIF (TM0.EQ.0.D0) THEN
          TEV(I) = TIME + STEP(I1)
          ITERM = 1
*
*       Increase TEV for standard case. 
      ELSE
          DT = MIN(TGR,TBR,TM0)
          IF(RTMS.GT.RL1)THEN
*       Estimate time of Roche-lobe filling for MS star. 
             DT1 = DTMS*(RL1 - RADIUS(J1))/(RTMS - RADIUS(J1))
             DT = MIN(DT1,DT)
          ENDIF
          TEV0(I) = TEV(I)
          TEV(I) = TEV(I) + DT
          TEV(I1) = TEV(I) + STEP(I1)
          TEV(I2) = TEV(I1)
      END IF
*
   30 RETURN
*
      END
