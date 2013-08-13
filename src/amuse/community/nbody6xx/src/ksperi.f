      SUBROUTINE KSPERI(IPAIR)
*
*
*       Pericentre KS variables.
*       ------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
*
*
*       Save current time and initialize integration indicator.
      TIME0 = TIME
      ITIME = 0
      I1 = 2*IPAIR - 1
      ICM = N + IPAIR
*
*       Re-calculate TDOT2 which may have misleading value.
    1 TD2 = 0.0D0
      DO 5 K = 1,4
          TD2 = TD2 + U0(K,IPAIR)*UDOT(K,IPAIR)
    5 CONTINUE
      TD2 = 2.0*TD2
*
*       Define auxiliary quantities and obtain the eccentricity.
      SEMI = -0.5D0*BODY(ICM)/H(IPAIR)
      ZETA = 1.0 - R(IPAIR)/SEMI
      PSI = TD2/SQRT(BODY(ICM))
      ECC = SQRT(ZETA**2 + PSI**2/SEMI)
*
*       Avoid nearly circular orbits (undefined pericentre).
      IF (ECC.LT.0.0001) GO TO 30
*
*       Distinguish between near-collision, elliptic & hyperbolic case.
      Q = (PSI/ECC)**2
      IF (ZETA.GT.0.0D0.AND.ABS(Q/SEMI).LT.0.1) THEN
*       Expansion of THETA/(SIN(THETA)*COS(THETA)) - 1 for Q = SIN(THETA)**2.
          SN1 = 1.0
          S = 0.0D0
          Q = Q/SEMI
          DO 8 IT = 1,10
              SN1 = Q*SN1*FLOAT(2*IT)/FLOAT(2*IT + 1)
              S = S + SN1
    8     CONTINUE
          S = S + SN1*Q*0.9D0/(1.0 - Q)
          DT = (R(IPAIR)*ZETA - PSI**2 + ZETA*S*SEMI)*PSI/
     &                                          (ECC**2*SQRT(BODY(ICM)))
      ELSE IF (SEMI.GT.0.0D0) THEN
*       Determine the eccentric anomaly with respect to pericentre (0,PI).
          THETA = ATAN2(ABS(PSI)/SQRT(SEMI),ZETA)
*       Obtain pericentre time interval from Kepler's equation.
          DT = SEMI*SQRT(SEMI/BODY(ICM))*(THETA - ABS(PSI)/SQRT(SEMI))
      ELSE IF (SEMI.LT.0.0D0) THEN
          A1 = PSI/(ECC*SQRT(ABS(SEMI)))
*       Use EXP(F) = SINH(F) + COSH(F) to obtain the eccentric anomaly THETA.
          A2 = ABS(A1) + SQRT(A1**2 + 1.0D0)
          THETA = LOG(A2)
          IF (A1.LT.0.0D0) THETA = -THETA
          A0 = ABS(SEMI)
          DT = A0*SQRT(A0/BODY(ICM))*(ABS(PSI)/SQRT(A0) - THETA)
      END IF
*
*       Re-define current time (NB! not quantized; T0(I1) may be << TIME).
      TIME = TIME - DT
*
*       Integrate backwards for perturbed motion (reflection gives errors).
      IF (LIST(1,I1).GT.0.AND.ITIME.EQ.0.AND.DT.GT.STEP(I1)) THEN
          TIME = TIME0
          IMOD = KSLOW(IPAIR)
          ZMOD = FLOAT(ISLOW(IMOD))
          IPH = IPHASE
          IPHASE = -1
*
*       Integrate step by step if interval is too large (note IPHASE < 0).
   10     IF (DT.GT.STEP(I1)) THEN
              TIME = TIME - STEP(I1)
              DT = DT - STEP(I1)
              H0(IPAIR) = H(IPAIR)
              Z = -0.5D0*H(IPAIR)*DTAU(IPAIR)**2
              CALL STUMPF(IPAIR,Z)
              DTAU(IPAIR) = -ABS(DTAU(IPAIR))
              CALL KSINT(I1)
              DTU = DTAU(IPAIR)
*       Use negative DTU and treat STEP as positive (not used elsewhere).
              STEP(I1) = ((ONE6*TDOT3(IPAIR)*DTU + 0.5*TDOT2(IPAIR))*DTU
     &                                                   + R(IPAIR))*DTU
              STEP(I1) = -ZMOD*STEP(I1)
              ITIME = ITIME + 1
              IF (ITIME.LT.200) GO TO 10
          END IF
*
          ITIME = ITIME + 1
          DTU = DT/(R(IPAIR)*ZMOD)
          DTU0 = DTAU(IPAIR)
          ITER = 0
*       Determine the regularized step by Newton-Raphson iteration (DT > 0).
   20     Y0 = DT - ZMOD*((ONE6*TDOT3(IPAIR)*DTU +
     &                             0.5*TDOT2(IPAIR))*DTU + R(IPAIR))*DTU
          YPR = -((0.5*TDOT3(IPAIR)*DTU + TDOT2(IPAIR))*DTU + R(IPAIR))
          YPR = ZMOD*YPR
          DTU = DTU - Y0/YPR
          DT1 = ((ONE6*TDOT3(IPAIR)*DTU + 0.5*TDOT2(IPAIR))*DTU +
     &                                                     R(IPAIR))*DTU
          DT1 = ZMOD*DT1
          ITER = ITER + 1
          IF (ABS(DT - DT1).GT.1.0D-10*STEP(I1).AND.ITER.LT.5) GO TO 20
*
*       Integrate back to pericentre using temporary indicator < 0 for exit.
          TIME = TIME - DT
          DTAU(IPAIR) = -DTU
*       Re-initialize Stumpff functions.
          H0(IPAIR) = H(IPAIR)
          Z = -0.5D0*H(IPAIR)*DTAU(IPAIR)**2
          CALL STUMPF(IPAIR,Z)
          CALL KSINT(I1)
*       Update energy and set positive step in case no KS initialization.
          H0(IPAIR) = H(IPAIR)
          DTAU(IPAIR) = ABS(DTU0)
          STEP(I1) = ZMOD*DT
          IPHASE = IPH
*       Use reflection procedure to improve provisional pericentre.
          GO TO 1
*       Note: typically 2 iterations and final TDOT2 > 0.
      END IF
*
*       Specify transformation coefficients (Seppo Mikkola's procedure).
      IF (ZETA.GE.0.0) THEN
          XC = SQRT(0.5D0 + 0.5D0*ZETA/ECC)
          YS = PSI/(ECC*XC*SQRT(BODY(ICM)))
      ELSE
*       Employ well behaved expressions for R > A (SM 29/5/97).
          XC = 0.5*ABS(PSI)/(SQRT(SEMI)*ECC)/SQRT(0.5D0-0.5D0*ZETA/ECC)
*       Avoid division by small XC near apocentre (ZETA < 0 only).
          YS = 2.0*SQRT(SEMI/BODY(ICM)*(0.5D0 - 0.5D0*ZETA/ECC))
          IF (PSI.LT.0.0) YS = -YS
      END IF
*
      ZZ = BODY(ICM)/(4.0*SEMI)
      R(IPAIR) = 0.0D0
      TDOT2(IPAIR) = 0.0D0
*
*       Generate analytical solutions for U & UDOT using old U0 & UDOT.
      T0(I1) = TIME
      DO 25 K = 1,4
          U(K,IPAIR) = U0(K,IPAIR)*XC - UDOT(K,IPAIR)*YS
          UDOT(K,IPAIR) = U0(K,IPAIR)*YS*ZZ + UDOT(K,IPAIR)*XC
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0D0*U(K,IPAIR)*UDOT(K,IPAIR)
   25 CONTINUE
*
*       Do not allow R' < 0 for repeated pericentre in routine KSINT.
      IF (TDOT2(IPAIR).LT.0.0D0) THEN
          TDOT2(IPAIR) = 0.0D0
      END IF
*
*       Predict c.m. coordinates & velocities.
      IF (ABS(DT).LT.STEP(ICM).AND.TIME.GT.0.0) THEN
          CALL XVPRED(ICM,0)
      END IF
*
   30 RETURN
*
      END
