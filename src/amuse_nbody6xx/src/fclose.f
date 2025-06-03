      SUBROUTINE FCLOSE(I,NNB)
*
*
*       Force & first derivative from close bodies.
*       -------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  A(9)
*
*
*       Initialize F & FDOT for body #I.
      DO 10 K = 1,3
          F(K,I) = 0.0D0
          FDOT(K,I) = 0.0D0
   10 CONTINUE
*
*       Obtain F & FDOT due to NNB members of JLIST.
      DO 50 L = 1,NNB
          J = JLIST(L)
          IF (J.EQ.I) GO TO 50
*
          RIJ2 = 0.0D0
          RIJDOT = 0.0D0
          DO 20 K = 1,3
              A(K) = X(K,J) - X(K,I)
              A(K+3) = XDOT(K,J) - XDOT(K,I)
              RIJ2 = RIJ2 + A(K)**2
              RIJDOT = RIJDOT + A(K)*A(K+3)
   20     CONTINUE
          A(8) = BODY(J)/(RIJ2*SQRT(RIJ2))
          A(9) = 3.0*RIJDOT/RIJ2
*
*       Set approximate F & FDOT to be used by body #I in FPOLY2.
          DO 30 K = 1,3
              F(K,I) = F(K,I) + A(K)*A(8)
              FDOT(K,I) = FDOT(K,I) + (A(K+3) - A(K)*A(9))*A(8)
   30     CONTINUE
   50 CONTINUE
*
*       Initialize time of last force (prevents prediction in FPOLY2).
      T0(I) = TIME
*
      RETURN
*
      END
