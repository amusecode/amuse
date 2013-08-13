      SUBROUTINE TPERI(SEMI,Q,UPR,MB,DT)
*
*
*       Pericentre time for KS motion.
*       ------------------------------
*
*       Routine for calculating the pericentre time
*       from the KS coordinates Q and derivatives UPR.
*       SEMI = semi-major axis of dominant motion.
*       Q = KS coordinates of relative position vector.
*       UPR = derivatives of Q in standard KS formulation.
*       MB = combined mass of the two bodies.
*       DT = pericentre time interval (DT < 0 before peri).
*       Algorithm due to Seppo Mikkola (1991).
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      REAL*8  Q(4),UPR(4)
*
*
*       Set scalar KS distance, velocity and radial velocity.
      R = Q(1)**2 + Q(2)**2 + Q(3)**2 + Q(4)**2
      V2 = UPR(1)**2 + UPR(2)**2 + UPR(3)**2 + UPR(4)**2
      VR = 2.0D0*(Q(1)*UPR(1) + Q(2)*UPR(2) + Q(3)*UPR(3) + Q(4)*UPR(4))
*
*       Define auxiliary variables.
      PSI = VR/SQRT(MB)
      ZETA = 4.0D0*V2/MB - 1.0D0
*
*       Form semi-major axis & eccentricity.
*     ALPHA = 2.0D0*(1.0D0 - 2.0D0*V2/MB)/R
*     SEMI = 1.0D0/ALPHA
*
*       Employ regular value of semi-major axis determined by routine EREL.
      ALPHA = 1.0D0/SEMI
      ECC = SQRT(ZETA*ZETA + ALPHA*PSI*PSI)
*
*       Distinguish between near-collision, elliptic & hyperbolic case.
      Q1 = (PSI/ECC)**2
      IF (ZETA.GT.0.0D0.AND.ABS(Q1/SEMI).LT.0.1) THEN
*       Expansion of THETA/(SIN(THETA)*COS(THETA)) - 1 for Q = SIN(THETA)**2.
          SN1 = 1.0
          S = 0.0D0
          Q1 = Q1/SEMI
          DO 1 IT = 1,10
              SN1 = Q1*SN1*FLOAT(2*IT)/FLOAT(2*IT + 1)
              S = S + SN1
    1     CONTINUE
          S = S + SN1*Q1*0.9D0/(1.0 - Q1)
          DT = (R*ZETA - PSI**2 + ZETA*S*SEMI)*PSI/(ECC**2*SQRT(MB))
      ELSE IF (SEMI.GT.0.0D0) THEN
*       Determine the eccentric anomaly with respect to pericentre (-PI,PI).
          THETA = DATAN2(PSI/SQRT(SEMI),ZETA)
*       Obtain pericentre time interval from Kepler's equation.
          DT = SEMI*SQRT(SEMI/MB)*(THETA - PSI/SQRT(SEMI))
      ELSE IF (SEMI.LT.0.0D0) THEN
*       Hyperbolic case.
          A1 = PSI/(ECC*SQRT(ABS(SEMI)))
*       Use EXP(F) = SINH(F) + COSH(F) to obtain the eccentric anomaly THETA.
          A2 = ABS(A1) + SQRT(A1**2 + 1.0D0)
          THETA = LOG(A2)
          IF (A1.LT.0.0D0) THETA = -THETA
          A0 = ABS(SEMI)
          DT = A0*SQRT(A0/MB)*(PSI/SQRT(A0) - THETA)
      END IF
*
      RETURN
*
      END
