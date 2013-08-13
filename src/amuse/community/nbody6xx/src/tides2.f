      SUBROUTINE TIDES2(QPERI,M1,M2,S1,S2,KSTAR,ECC,KG,WS,QS,DE2,DE3)
*
*
*       Tidal energy loss for interacting bodies.
*       -----------------------------------------
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  SIZE(2),RATIO(2),ZETA(2),PSI2(2),PSI3(2),DE2(2),DE3(2),
     &        EFAC2,EFAC3,WS(2),QS(2),M1,M2,BODY(2)
      INTEGER KSTAR(2),IP(2),KG(2)
*
*
*       Copy masses & radii to local arrays.
      BODY(1) = M1
      BODY(2) = M2
      SIZE(1) = S1
      SIZE(2) = S2
      EFAC = 2.0/(1.0 + ECC)
*
*       Assign the appropriate polytropic reference index (n = 3, 2 or 3/2).
      DO 5 K = 1,2
          II = 3
          IF (KSTAR(K).EQ.0) II = 1
          IF (KSTAR(K).GE.10) II = 2
          IP(K) = II
    5 CONTINUE
*
*       Form dimensionless ratios & tidal efficiency factors for each star.
      DO 10 K = 1,2
          L = 3 - K
          RATIO(K) = SIZE(K)/QPERI
          IF (RATIO(K).GT.0.1) THEN
              ICASE = IP(K)
              ZETA(K) = SQRT((BODY(K)/(BODY(K) + BODY(L)))/RATIO(K)**3)
              IF (ZETA(K).GT.2.0) THEN
                  ALPHA = 0.5 + 0.25*(0.5*(ZETA(K) - 2.0))**1.5
              ELSE
                  ALPHA = 0.5
              END IF
              ZETA(K) = ZETA(K)*EFAC**ALPHA
              IF (KG(K).EQ.0) THEN
                  PSI2(K) = EFAC2(ZETA(K),ICASE)
                  PSI3(K) = EFAC3(ZETA(K),ICASE)
              ELSE
                  PSI2(K) = QS(1)*EFAC2(WS(1)*ZETA(K),ICASE)
                  PSI3(K) = QS(2)*EFAC3(WS(2)*ZETA(K),ICASE)
              END IF
          END IF
   10 CONTINUE
*
*       Obtain energy loss due to second & third order modes.
      DO 20 K = 1,2
          L = 3 - K
*         IF (RATIO(K).GT.0.1) THEN
          IF (RATIO(K).GT.0.5) THEN    ! awaiting bug fix (too large values)
              R2 = RATIO(K)**2
              R6 = R2*R2*R2*BODY(L)**2/SIZE(K)
              DE2(K) = R6*PSI2(K)
              DE3(K) = R6*R2*PSI3(K)
          ELSE
              DE2(K) = 0.0
              DE3(K) = 0.0
          END IF
   20 CONTINUE
*
      RETURN
*
      END
