      SUBROUTINE NEWSYS(X,XD,M,NP,ENERGY,GAM)
*
*
*       Total energy of subsystem.
*       --------------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  M(*),X(*),XD(*)
*
*
*       Calculate total energy in physical units.
      T = 0.0D0
      V = 0.0D0
      DO 10 I = 1,NP
          K1 = 3*(I - 1) + 1
          K2 = K1 + 1
          K3 = K2 + 1
          T = T + 0.5D0*M(I)*(XD(K1)**2 + XD(K2)**2 + XD(K3)**2)
          IF (I.EQ.NP) GO TO 10
          J1 = I + 1
          DO 9 J = J1,NP
              KI = 3*(I - 1)
              KJ = 3*(J - 1)
              R2 = 0.0D0
              DO 8 K = 1,3
                  KI = KI + 1
                  KJ = KJ + 1
                  R2 = R2 + (X(KI) - X(KJ))**2
    8         CONTINUE
              V = V - M(I)*M(J)/SQRT(R2)
    9     CONTINUE
   10 CONTINUE
*
      ENERGY = T + V
      GAM = T - V
*
      RETURN
*
      END
