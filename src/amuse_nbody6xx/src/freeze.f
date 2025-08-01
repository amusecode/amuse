      SUBROUTINE FREEZE(IPAIR)
*
*
*       Partial reflection of KS binary.
*       --------------------------------
*
      INCLUDE 'common6.h'
      SAVE NKSREF
      DATA NKSREF /0/
*
*
*       Set c.m. & component indices and form semi-major axis.
      ICM = N + IPAIR
      I2 = 2*IPAIR
      I1 = I2 - 1
      SEMI = -0.5D0*BODY(ICM)/H(IPAIR)
*
*       Define auxiliary quantities and obtain the eccentricity.
      ZETA = 1.0 - R(IPAIR)/SEMI
      PSI = TDOT2(IPAIR)/SQRT(BODY(ICM))
      ECC = SQRT(ZETA**2 + PSI**2/SEMI)
*
*       Skip reflection of hyperbolic orbit (just in case).
      IF (ECC.GE.1.0) GO TO 100
*
*       Determine the eccentric anomaly with respect to pericentre (0,PI).
      THETA = ATAN2(ABS(PSI)/SQRT(SEMI),ZETA)
*
*       Obtain total reflection time from Kepler's equation.
      DT = 2.0D0*SEMI*SQRT(SEMI/BODY(ICM))*(THETA - ABS(PSI)/SQRT(SEMI))
*
*       Specify regularized time (based on Baumgarte & Stielel, 1974).
      DTU = -2.0D0*(H(IPAIR)*DT + TDOT2(IPAIR))/BODY(ICM)
*
*       Skip reflection near pericentre.
      IF (DTU.LT.4.0*DTAU(IPAIR)) GO TO 100
*
*       Predict c.m. coordinates and resolve unreflected KS components.
      CALL XVPRED(ICM,0)
*
*       Specify transformation coefficients (Mikkola's proceure).
      XC = ZETA/ECC
      YS = 2.0D0*PSI/(ECC*SQRT(BODY(ICM)))
      ZZ = BODY(ICM)/(4.0*SEMI)
      R(IPAIR) = 0.0D0
*
*       Generate analytical solutions for U & UDOT using old U & UDOT.
      DO 10 K = 1,4
          U(K,IPAIR) = U0(K,IPAIR)*XC - UDOT(K,IPAIR)*YS
          UDOT(K,IPAIR) = U0(K,IPAIR)*YS*ZZ + UDOT(K,IPAIR)*XC
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
*       Consistent R = U*U for stabilization procedure.
   10 CONTINUE
*
*       Copy significant perturbers for tidal correction (R-dependent).
      CORR = 2.0 + 2.0*(SEMI - R(IPAIR))/(SEMI*(1.0 + ECC))
      RCRIT2 = CMSEP2*(CORR*R(IPAIR))**2
      NNB1 = LIST(1,I1) + 1
      NP = 0
      DO 20 L = 2,NNB1
          J = LIST(L,I1)
          RIJ2 = (X(1,ICM) - X(1,J))**2 + (X(2,ICM) - X(2,J))**2 +
     &                                    (X(3,ICM) - X(3,J))**2
          IF (RIJ2.LT.RCRIT2) THEN
              NP = NP + 1
              JPERT(NP) = J
          END IF
   20 CONTINUE
*
*       Ensure that at least one perturber is retained.
      IF (NP.EQ.0) THEN
          NP = 1
          JPERT(NP) = LIST(2,I1)
      END IF
*
*       Reverse second time derivative and specify unperturbed motion.
      TDOT2(IPAIR) = -TDOT2(IPAIR)
      LIST(1,I1) = 0
*
*       Set new step in physical units and increase counter.
      STEP(I1) = DT
      NKSREF = NKSREF + 1
*
*       Check minimum two-body distance.
      DMIN2 = MIN(DMIN2,SEMI*(1.0D0 - ECC))
*
*       Obtain potential energy with respect to the components.
      JLIST(1) = I1
      JLIST(2) = I2
      CALL NBPOT(2,NP,POT1)
*
*       Transform to reflected coordinates and repeat the summation.
      CALL KSRES(IPAIR,J1,J2,RIJ2)
      CALL NBPOT(2,NP,POT2)
*
*       Form correction factor due to the tidal potential.
      VI2 = X0DOT(1,ICM)**2 + X0DOT(2,ICM)**2 + X0DOT(3,ICM)**2
      CORR = 1.0 + 2.0*(POT2 - POT1)/(BODY(ICM)*VI2)
      ETCORR = ETCORR + (POT2 - POT1)
      IF (CORR.GT.0.0D0) THEN
          CORR = SQRT(CORR)
      ELSE
          CORR = 0.0D0
      END IF
*
*       Modify the c.m. velocity.
      DO 30 K = 1,3
          X0DOT(K,ICM) = CORR*X0DOT(K,ICM)
   30 CONTINUE
*
  100 RETURN
*
      END
