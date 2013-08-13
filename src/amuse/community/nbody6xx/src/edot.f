      SUBROUTINE EDOT(IPAIR,J,SEMI,ECC,ECCDOT)
*
*
*       Eccentricity derivative due to perturber.
*       -----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  XREL(3),VREL(3),FP(3)
*
*
*       Set relative coordinates & velocities and initialize perturbation.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      DO 1 K = 1,3
          XREL(K) = X(K,I1) - X(K,I2)
          VREL(K) = XDOT(K,I1) - XDOT(K,I2)
          FP(K) = 0.0D0
    1 CONTINUE
*
*       Obtain perturbation on KS components due to body #J.
      L = I1
    5 A1 = X(1,J) - X(1,L)
      A2 = X(2,J) - X(2,L)
      A3 = X(3,J) - X(3,L)
      RIJ2 = A1**2 + A2**2 + A3**2
      A4 = BODY(J)/(RIJ2*SQRT(RIJ2))
      IF (L.EQ.I2) A4 = -A4
      FP(1) = FP(1) + A1*A4
      FP(2) = FP(2) + A2*A4
      FP(3) = FP(3) + A3*A4
      IF (L.EQ.I1) THEN
          L = I2
          GO TO 5
      END IF
*
*       Form scalar products R*V, R*F and V*F.
      RV = 0.0
      RF = 0.0
      VF = 0.0
      DO 10 K = 1,3
          RV = RV + XREL(K)*VREL(K)
          RF = RF + XREL(K)*FP(K)
          VF = VF + VREL(K)*FP(K)
   10 CONTINUE
*
*       Evaluate time derivative of eccentricity (Douglas Heggie 30/8/96).
      ECCDOT = (SEMI**2*(1.0 - ECC**2) - R(IPAIR)**2)*VF + RV*RF
      ECCDOT = ECCDOT/((BODY(I1) + BODY(I2))*SEMI*ECC)
*
*     WRITE (6,15)  IPAIR, KSTAR(N+IPAIR), ECC, ECCDOT, RV, RF, VF
*  15 FORMAT (' EDOT:    KS K* E ED RV RF VF ',2I4,F7.3,1P,4E10.2)
*
      RETURN
*
      END
