      SUBROUTINE CHPOT(DP)
*
*
*       Potential energy correction due to chain.
*       -----------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
*
*
*       Obtain current global coordinates of chain components.
      CALL XCPRED(0)
*
*       Consider contributions from all active perturbers.
      DP = 0.0D0
      JDUM = 0
      NNBC = LISTC(1) + 1
*
      DO 10 L = 2,NNBC
*       Subtract potential energy due to chain c.m.
          J = LISTC(L)
*       Replace any regularized c.m. body by individual components.
          IF (J.GT.N) THEN
              JDUM = 2*(J - N) - 1
              J = JDUM
          END IF
    2     A1 = X(1,ICH) - X(1,J)
          A2 = X(2,ICH) - X(2,J)
          A3 = X(3,ICH) - X(3,J)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          DP = DP - BODY(ICH)*BODY(J)/SQRT(RIJ2)
*
*       Add individual interactions to obtain differential correction.
          DO 5 K = 1,NCH
              A1 = XC(1,K) - X(1,J)
              A2 = XC(2,K) - X(2,J)
              A3 = XC(3,K) - X(3,J)
              RIJ2 = A1*A1 + A2*A2 + A3*A3
              DP = DP + BODYC(K)*BODY(J)/SQRT(RIJ2)
    5     CONTINUE
*
*       Check for possible second KS component and restore dummy index.
          IF (J.EQ.JDUM) THEN
              J = JDUM + 1
              JDUM = 0
              GO TO 2
          END IF
   10 CONTINUE
*
      RETURN
*
      END
