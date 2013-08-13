      SUBROUTINE TIDES(QPERI,M1,M2,S1,S2,KSTAR,DE)
*
*
*       Tidal energy loss for interacting stars.
*       ----------------------------------------
*
      REAL*8  SIZE(2),RATIO(2),ZETA(2),PSI2(2),PSI3(2),DE(2),
     &        QPERI,M1,M2,EFAC2,EFAC3,S1,S2,BODY(2)
      INTEGER KSTAR(2),IP(2)
*
*
*       Copy masses & radii to local arrays.
      BODY(1) = M1
      BODY(2) = M2
      SIZE(1) = S1
      SIZE(2) = S2
*
*       Assign the appropriate polytropic reference index (NP = 3/2 or 3).
      DO 5 K = 1,2
          KST = KSTAR(K)
          IF (KST.EQ.3.OR.KST.EQ.5.OR.KST.EQ.6) THEN
              IP(K) = 1
          ELSE
              IP(K) = 3
          END IF
    5 CONTINUE
*
*       Form dimensionless ratios & tidal efficiency factors for each star.
      DO 10 K = 1,2
          L = 3 - K
          RATIO(K) = SIZE(K)/QPERI
          IF (RATIO(K).GT.0.1) THEN
              ICASE = IP(K)
              ZETA(K) = SQRT((BODY(K)/(BODY(K) + BODY(L)))/RATIO(K)**3)
              PSI2(K) = EFAC2(ZETA(K),ICASE)
              PSI3(K) = EFAC3(ZETA(K),ICASE)
          END IF
   10 CONTINUE
*
*       Obtain energy loss due to second & third order modes.
      DO 20 K = 1,2
          L = 3 - K
          IF (RATIO(K).GT.0.1) THEN
              R2 = RATIO(K)**2
              R6 = R2*R2*R2
              DE(K) = R6*(PSI2(K) + R2*PSI3(K))*BODY(L)**2/SIZE(K)
          ELSE
              DE(K) = 0.0
          END IF
   20 CONTINUE
*
      RETURN
*
      END
