      SUBROUTINE SETUP
*
*
*       Generation of initial coordinates & velocities.
*       -----------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  A(8)
      REAL*8  RAN2
*
*
*       Choose between uniform density and Plummer model.
      KDUM = IDUM1
      IF (KZ(5).GT.0) GO TO 20
*
*       Set up a uniform spherical system.
      DO 10 I = 1,N
    1     A(1) = 0.0D0
          DO 2 K = 1,3
              A(K+1) = 2.0*RAN2(KDUM) - 1.0
              A(1) = A(1) + A(K+1)**2
    2     CONTINUE
          IF (A(1).GT.1.0) GO TO 1
*
    4     A(5) = 0.0D0
          DO 5 K = 1,3
              A(K+5) = 2.0*RAN2(KDUM) - 1.0
              A(5) = A(5) + A(K+5)**2
    5     CONTINUE
          IF (A(5).GT.1.0) GO TO 4
*
          DO 8 K = 1,3
*             X(K,I) = A(1)*A(K+1)
*       Density proportional to 1/R**2.
              X(K,I) = A(K+1)
*       Constant density.
              XDOT(K,I) = A(K+5)
*       Isotropic velocities (magnitude randomized; no radial dependence).
    8     CONTINUE
   10 CONTINUE
*
      GO TO 60
*
*       Initialize centre of mass terms.
   20 DO 25 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   25 CONTINUE
*
*       Generate initial conditions from Plummer model (A & A 37, 183).
      DO 40 I = 1,N
   30     A(1) = RAN2(KDUM)
          IF (A(1).LT.1.0D-10) GO TO 30
          RI = (A(1)**(-0.6666667) - 1.0)**(-0.5)
*       Reject distant particles.
          IF (RI.GT.10.0) GO TO 30
*
          A(2) = RAN2(KDUM)
          A(3) = RAN2(KDUM)
          X(3,I) = (1.0 - 2.0*A(2))*RI
          X(1,I) = SQRT(RI**2 - X(3,I)**2)*COS(TWOPI*A(3))
          X(2,I) = SQRT(RI**2 - X(3,I)**2)*SIN(TWOPI*A(3))
   32     A(4) = RAN2(KDUM)
          A(5) = RAN2(KDUM)
          A(6) = A(4)**2*(1.0 - A(4)**2)**3.5
          IF (0.1*A(5).GT.A(6)) GO TO 32
*
          A(8) = A(4)*SQRT(2.0)/(1.0 + RI**2)**0.25
          A(6) = RAN2(KDUM)
          A(7) = RAN2(KDUM)
          XDOT(3,I) = (1.0 - 2.0*A(6))*A(8)
          XDOT(1,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*COS(TWOPI*A(7))
          XDOT(2,I) = SQRT(A(8)**2 - XDOT(3,I)**2)*SIN(TWOPI*A(7))
*
*       Accumulate c.m. terms.
          DO 35 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   35     CONTINUE
   40 CONTINUE
*
*       Scale coordinates & velocities to analytical expectation values.
      SX = 1.5*TWOPI/16.0
      SV = SQRT(ZMASS/SX)
      DO 50 I = 1,N
          DO 45 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
              X(K,I) = SX*X(K,I)
              XDOT(K,I) = SV*XDOT(K,I)
   45     CONTINUE
   50 CONTINUE
*
*
*       Save random number sequence in COMMON for future use.
   60 IDUM1 = KDUM
*
      RETURN
*
      END
