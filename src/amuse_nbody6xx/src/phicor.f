      SUBROUTINE PHICOR(I,DPHI1,DPHI2)
*
*
*       Differential potential corrections.
*       -----------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Distinguish between single particle and (un)perturbed binary c.m.
      DPHI1 = 0.0
      DPHI2 = 0.0
      IF (I.LE.N) THEN
          NNB = LIST(1,I)
*       Search neighbour list backwards until the last c.m. particle.
          DO 10 L = NNB+1,2,-1
              J = LIST(L,I)
              IF (J.LE.N) GO TO 40
              RIJ2 = 0.0
              DO 5 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
    5         CONTINUE
*       Perform differential correction inside standard c.m. separation.
              IF (RIJ2.LT.CMSEP2*RMIN2) THEN
                  RJ1 = 0.0
                  RJ2 = 0.0
                  J1 = 2*(J - N) - 1
                  J2 = J1 + 1
                  DO 8 K = 1,3
                      RJ1 = RJ1 + (X(K,I) - X(K,J1))**2
                      RJ2 = RJ2 + (X(K,I) - X(K,J2))**2
    8             CONTINUE
                  DPHI1 = DPHI1 - BODY(J)/SQRT(RIJ2) +
     &                           BODY(J1)/SQRT(RJ1) + BODY(J2)/SQRT(RJ2)
              END IF
   10     CONTINUE
      ELSE
          I1 = 2*(I - N) - 1
          I2 = I1 + 1
          NP = LIST(1,I1)
*       Loop over all perturbers (only c.m. approximation for singles).
          DO 30 L = 2,NP+1
              J = LIST(L,I1)
              RIJ2 = 0.0
              DO 12 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   12         CONTINUE
*       Treat case of close perturbing binary more carefully below.
              IF (J.LE.N.OR.RIJ2.GT.CMSEP2*RMIN2) THEN
                  RI1 = 0.0
                  RI2 = 0.0
                  DO 15 K = 1,3
                      RI1 = RI1 + (X(K,I1) - X(K,J))**2
                      RI2 = RI2 + (X(K,I2) - X(K,J))**2
   15             CONTINUE
*       Subtract c.m. interaction for each component (copied from GPU).
                  DPHI1 = DPHI1 - BODY(J)/SQRT(RIJ2) + BODY(J)/SQRT(RI1)
                  DPHI2 = DPHI2 - BODY(J)/SQRT(RIJ2) + BODY(J)/SQRT(RI2)
              ELSE
                  DPHI1 = DPHI1 - BODY(J)/SQRT(RIJ2)
                  DPHI2 = DPHI2 - BODY(J)/SQRT(RIJ2)
*       Include pair-wise interactions after subtracting the c.m.-c.m. terms.
                  J1 = 2*(J - N) - 1
                  J2 = J1 + 1
   25             RJ1 = 0.0
                  RJ2 = 0.0
                  DO 28 K = 1,3
                      RJ1 = RJ1 + (X(K,I1) - X(K,J1))**2
                      RJ2 = RJ2 + (X(K,I2) - X(K,J1))**2
   28             CONTINUE
*       Add each interaction separately (J1 followed by J2).
                  DPHI1 = DPHI1 + BODY(J1)/SQRT(RJ1)
                  DPHI2 = DPHI2 + BODY(J1)/SQRT(RJ2)
                  J1 = J1 + 1
                  IF (J1.EQ.J2) GO TO 25
              END IF
   30     CONTINUE
      END IF
*
   40 RETURN
*
      END
