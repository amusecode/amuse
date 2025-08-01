      SUBROUTINE BRAKE(IPAIR,DSEP,ECC1)
*
*
*       Orbital changes (GR, mass loss and/or tides).
*       ---------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  M1,M2
      SAVE NDIAG
      DATA NDIAG /0/
*
*
*       Set indices of KS components & c.m.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      I = N + IPAIR
      DSEP = DSEP/SU 
*
*       Set mass and radius of each star in solar units.
      M1 = BODY(I1)*ZMBAR
      M2 = BODY(I2)*ZMBAR
      R1 = RADIUS(I1)*SU
      R2 = RADIUS(I2)*SU
*
*       Obtain current semi-major axis & eccentricity and define new semi.
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
      SEMI1 = SEMI - DSEP 
*
*      Check for circularized case. 
      IF(KSTAR(I).LT.10.AND.ECC1.LT.0.002)THEN
         ECC1 = 0.001D0
         KSTAR(I) = 10
      ENDIF
*
*       Check collision/coalescence (IQCOLL=-2 skips ECOLL updating in COAL).
      RP = SEMI1*SU*(1.D0 - ECC1)
*
*       Include GR coalescence criterion for compact objects.
      KSX = MAX(KSTAR(I1),KSTAR(I2))
      IF (KSX.GE.13.AND.KZ(28).GT.0) THEN
         RCOAL = 6.0*BODY(I)/CLIGHT**2
      ELSE
         RCOAL = 100.0
      END IF
      IF (rank.eq.0.and.KZ(28).GT.1.AND.DSEP.GT.0.0) THEN
         WRITE (6,5)  TTOT, NAME(I1), KSX, ECC, SEMI*SU, DSEP/SEMI
    5    FORMAT (' BRAKE    T NAM KX* E A DA/A ',
     &                      F8.2,I6,I4,F9.5,1P,2E10.2)
      END IF
*
*       Check collision condition for degenerate objects.
      IF (RP.LE.(R1 + R2).OR.
     &   (KSX.GE.13.AND.R(IPAIR).LT.RCOAL)) THEN
         CALL KSPERI(IPAIR)
         IQCOLL = -2
         KSPAIR = IPAIR
         CALL CMBODY(R(IPAIR),2)
         GOTO 50
      ENDIF
      IF(DSEP.EQ.0.D0) GOTO 50
*
*       Include safety test on new semi-major axis.
      RCHCK = MIN(RADIUS(I1),RADIUS(I2)) 
      IF(SEMI1.LT.RCHCK) SEMI1 = RCHCK 
*
*       Transform to pericentre (R > A & unperturbed).
*     IF(R(IPAIR).GT.SEMI.AND.LIST(1,I1).EQ.0)THEN
*        CALL KSAPO(IPAIR)
*     ENDIF
*
*       Form square of regularized velocity.
      V20 = 0.0
      DO 10 K = 1,4
         V20 = V20 + UDOT(K,IPAIR)**2
   10 CONTINUE
*
*       Update binding energy and collision energy.
      HI = H(IPAIR)
      H(IPAIR) = -0.5*BODY(I)/SEMI1
      ZMU = BODY(I1)*BODY(I2)/BODY(I) 
      ECOLL = ECOLL + ZMU*(HI - H(IPAIR))
      EGRAV = EGRAV + ZMU*(HI - H(IPAIR))
*
*       Distinguish between update of eccentric binaries and standard case.
      IF(ABS(ECC-ECC1).GT.TINY)THEN
*       Change KS variables at original ecc and modify ECC at H = const.
         CALL EXPAND(IPAIR,SEMI)
         CALL KSPERI(IPAIR)
         CALL DEFORM(IPAIR,ECC,ECC1)
      ELSE
*       Specify KS coordinate & velocity scaling factors at arbitrary point.
         C2 = SQRT(SEMI1/SEMI)
*        V2 = 0.5*(BODY(I) + H(IPAIR)*SEMI1*(1.0 - ECC))
         V2 = 0.5*(BODY(I) + H(IPAIR)*R(IPAIR)*(SEMI1/SEMI))
         C1 = SQRT(V2/V20)
*
*       Re-scale KS variables to new energy with constant eccentricity.
         R(IPAIR) = 0.0D0
*        TDOT2(IPAIR) = 0.0D0
         DO 20 K = 1,4
            U(K,IPAIR) = C2*U(K,IPAIR)
            UDOT(K,IPAIR) = C1*UDOT(K,IPAIR)
            U0(K,IPAIR) = U(K,IPAIR)
            R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
*           TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0*U(K,IPAIR)*UDOT(K,IPAIR)
   20    CONTINUE
      END IF
*
*       Transform back to apocentre for standard unperturbed motion.
*     IF (LIST(1,I1).EQ.0) THEN
*        CALL KSAPO(IPAIR)
*     END IF
*
*       Include new initialization for perturbed orbit.
      IF (LIST(1,I1).GT.0) THEN
         IMOD = KSLOW(IPAIR)
         CALL KSPOLY(IPAIR,IMOD)
      ENDIF
*
*       Include some diagnostic output.
      IF(KSTAR(I).EQ.13)THEN
         NDIAG = NDIAG + 1
         IF(NDIAG.LT.100.OR.MOD(NDIAG,100).EQ.100)THEN
            RCOLL = RADIUS(I1) + RADIUS(I2)
            if(rank.eq.0)
     &      WRITE (6,25)  TTOT, IPAIR, M1, M2, R1, R2, R(IPAIR),
     &                    RCOLL
   25       FORMAT (' BRAKE    T KS M12 R12 R RCOLL ',
     &                         F10.4,I4,2F6.2,2F7.3,1P,2E10.2)
         ENDIF
      ENDIF
*
   50 RETURN
*
      END
