      SUBROUTINE KEPLER(I,DTI)
*
*
*       Step reduction of hierarchy.
*       ----------------------------
*
*
      INCLUDE 'common6.h'
*
*
*       Only examine hierarchical configurations (NP < 5 OR STEP < DTMIN).
      I1 = 2*(I - N) - 1
      IF (LIST(1,I1).GT.4.AND.DTI.GT.DTMIN) GO TO 30
*
*       Set list membership & Kepler period.
      NP1 = LIST(1,I1) + 1
      SEMI = -0.5D0*BODY(I)/H(I-N)
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
*
*       Consider the most dominant perturbers having comparable step to c.m.
      DO 20 L = 2,NP1
          J = LIST(L,I1)
          IF (STEP(J).GT.4.0*STEP(I)) GO TO 20 
          RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                  (X(3,I) - X(3,J))**2
          RIJ = SQRT(RIJ2)
          DT2 = 0.1*RIJ*SQRT(ETAI*RIJ/(BODY(I) + BODY(J)))
          DT = 0.25D0*SQRT(ETAI)*RIJ*TK/SEMI
*       Compare predicted c.m. step with conservative Kepler expressions.
          DT = MIN(DT,DT2)
          IF (DTI.LT.DT) GO TO 20 
          DTI = DT
*
*       Check whether to reduce step of dominant perturber.
          IF (STEP(J).LT.DT) GO TO 20
          IF (T0(J).EQ.TIME) THEN
              STEP(J) = 0.5D0*STEP(J)
              TIMENW(J) = T0(J) + STEP(J)
              GO TO 20
          END IF
          JCLOSE = J
   20 CONTINUE
*
   30 RETURN
*
      END
