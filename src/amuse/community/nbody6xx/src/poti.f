***
      SUBROUTINE POTI(I,POTJ)
*
*
*       Potential of one particle.
*       --------------------------
*
      INCLUDE 'common6.h'
*
*
*       Obtain the potential of body #I on host.
      POTJ = 0.D0
      KDUM = 0
      DO 30 JDUM = IFIRST,NTOT
          IF (JDUM.EQ.I) GO TO 30
          J = JDUM
          IF (J.GT.N) THEN
              JPAIR = J - N
*       Use c.m. approximation for unperturbed binary.
              IF (LIST(1,2*JPAIR-1).GT.0) THEN
                  KDUM = 2*JPAIR - 1
                  J = KDUM
              END IF
          END IF
   20     RIJ2 = 0.D0
          DO 25 K = 1,3
              RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   25     CONTINUE
          POTJ = POTJ - BODY(J)/SQRT(RIJ2)
          IF (J.EQ.KDUM) THEN
              J = J + 1
              GO TO 20
          END IF
   30 CONTINUE
*
      RETURN
      END
***
