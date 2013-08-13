      SUBROUTINE BRAKE3(IPAIR,ITERM)
*
*
*       Gravitational radiation of hierarchical BH/NS binary.
*       -----------------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/SLOW0/  RANGE,ISLOW(10)
      PARAMETER  (AURSUN=214.95D0)
      REAL*8  M1,M2,JORB
*
*
*       Define scalars and determine merger index.
      I = N + IPAIR
      I1 = 2*IPAIR - 1
      ITERM = 0
      IM = 0
      DO 1 K = 1,NMERGE
         IF (NAMEM(K).EQ.NAME(I)) IM = K
    1 CONTINUE
*
*       Find location of secondary (i.e. ghost) and form mass factors & SEMI.
      CALL FINDJ(I,I2,IM)
      IF (I2.LT.0) GO TO 30
      IF(KSTAR(I2).LT.10)THEN
         ITERM = 1
         GO TO 30
      ENDIF
      M1 = CM(1,IM)*ZMBAR
      M2 = CM(2,IM)*ZMBAR
      ZMB = CM(1,IM) + CM(2,IM)
      ZMU = CM(1,IM)*CM(2,IM)/ZMB
      SEMI = -0.5*ZMB/HM(IM)
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
*       Form inner eccentricity and angular momentum.
      SEP = SEMI*SU
      ECC2 = (1.0 - RB/SEMI)**2 + TD2**2/(ZMB*SEMI)
      ECC0 = SQRT(ECC2)
      OORB = TWOPI*SQRT((M1 + M2)/(SEP/AURSUN)**3)
      JORB = M1*M2/(M1 + M2)*SQRT(1.D0-ECC2)*SEP**2*OORB
*
*       Obtain angular momentum and eccentricity derivatives (SEMI < 10 R_s).
      CALL GRRAD(M1,M2,SEP,ECC0,JORB,DJGR,DELET)
*
*       Allow 2% angular momentum change (expressions adapted from MDOT).
      IF (ABS(DJGR).GT.0.D0) THEN
          DTGR = 0.02D0*JORB/ABS(DJGR)
          IF(DELET.GT.TINY.AND.ECC0.GT.0.0011D0)THEN
             DTGR = MIN(DTGR,0.05D0*ECC0/DELET)
          ENDIF
          DTGR = MAX(DTGR,100.D0)
*       Note that time interval is in units of yr (cf. grrad.f).
          DJORB = DJGR*DTGR
*       Terminate on large change of JORB which occurs at GR coalescence.
          IF (DJORB.GT.0.1*JORB) THEN
              ITERM = 1
              GO TO 30
          END IF
          JORB = JORB - DJORB
          JORB = MAX(JORB,1.D0)
          ECC = MAX(ECC0 - DELET*DTGR,0.001D0)
*       Convert to N-body units and update both TEV to same time (NB!).
          DTGR = DTGR/(1.0D+06*TSTAR)
          TEV0(I1) = TIME
          TEV0(I2) = TIME
          TEV(I1) = TIME + DTGR
          TEV(I2) = TIME + DTGR
      ELSE
          DT = MIN(TEV(I1)-TEV0(I1),TEV(I2)-TEV0(I2))
          DT = MAX(DT,0.1D0)
          TEV0(I1) = TIME
          TEV0(I2) = TIME
          TEV(I1) = TIME + DT
          TEV(I2) = TIME + DT
          GO TO 30
      END IF
*
*       Obtain new semi-major axis (in N-body units) from angular momentum.
      SEMI1 = (M1 + M2)*JORB*JORB/
     &                  ((M1*M2*TWOPI)**2*AURSUN**3*(1.D0-ECC**2))
      SEMI1 = SEMI1/SU
*
*       Update inner binding energy.
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
*       Specify new peri- or apo-centre at turning-point.
      IF (RB.LT.SEMI1) THEN
          RNEW = SEMI1*(1.0 - ECC)
          EFAC = (1.0 - ECC)/(1.0 - ECC0)
      ELSE
          RNEW = SEMI1*(1.0 + ECC)
          EFAC = (1.0 + ECC)/(1.0 + ECC0)
      END IF
      C1 = SQRT(EFAC)
      V2 = 0.5*(ZMB + HM(IM)*RNEW)
      IF (V2.LE.0.D0) THEN
         C2 = 1.0D-06
      ELSE
         C2 = SQRT(V2/V20)
      END IF
*
*       Re-scale KS variables at constant energy to new eccentricity (ECC).
      RB = 0.0D0
      DO 20 K = 1,4
          UM(K,IM) = C1*UM(K,IM)
          UMDOT(K,IM) = C2*UMDOT(K,IM)
          RB = RB + UM(K,IM)**2
   20 CONTINUE
*
*       Rectify the hierarchical KS variables.
      CALL HIRECT(IM)
*
*       Include GR coalescence criterion for compact objects.
      KSX = MAX(KSTAR(I1),KSTAR(I2))
      IF (KSX.GE.13.AND.KZ(28).GT.0) THEN
         RCOAL = 6.0*ZMB/CLIGHT**2
      ELSE
         RCOAL = 100.0
      END IF
*
      WRITE (6,25)  NAME(I1), TOFF+TIME, ECC, SEMI1, DTGR, RCOAL
   25 FORMAT (' BRAKE3    NM T E A DTGR RCOAL',I8,F9.2,F8.4,1P,3E10.2)
*
*       Check termination criterion for coalescence.
      IF (RB.LT.RCOAL) THEN
          ITERM = 1
      END IF
   30 CONTINUE
*
      RETURN
*
      END
