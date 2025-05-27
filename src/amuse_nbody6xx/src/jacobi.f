      SUBROUTINE JACOBI(NESC)
*
*
*       Jacobi escape criterion.
*       ------------------------
*
      INCLUDE 'common6.h'
*
*
*       Specify escape energy (tidal field or isolated system).
      IF (KZ(14).GT.0.AND.KZ(14).LT.2) THEN
          ECRIT = -1.5*(TIDAL(1)*ZMASS**2)**0.333
      ELSE
          ECRIT = 0.0
      END IF
*       Define current mean square velocity (= 0.5 initially).
      V2M = 0.5*ZMASS/RSCALE
*
*       Count all escapers.
      NESC = 0
      DO 60 I = IFIRST,NTOT
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                                   (X(3,I) - RDENS(3))**2
          VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
          IF (RI2.GT.RSCALE**2) THEN
              VC2 = ZMASS/SQRT(RI2) + ECRIT
              IF (VI2.LT.VC2) GO TO 60
          ELSE
*       Skip velocities below sqrt(2) times equilibrium value.
              IF (VI2.LT.2.0*V2M + ECRIT) GO TO 60
          END IF
          POTI = 0.0
          DO 55 J = IFIRST,NTOT
              IF (J.EQ.I) GO TO 55
              RIJ2 = 0.0
              DO 52 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   52         CONTINUE
              POTI = POTI + BODY(J)/SQRT(RIJ2)
   55     CONTINUE
          EI = 0.5*VI2 - POTI
          IF (EI.GT.ECRIT) NESC = NESC + 1
   60 CONTINUE
*
      RETURN
*
      END
