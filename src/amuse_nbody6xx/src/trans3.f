      SUBROUTINE TRANS3(KDUM)
*
*
*       Transformation of physical & KS variables.
*       ------------------------------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      COMMON/AZREG/  TIME3,TMAX,Q(8),P(8),R1,R2,R3,ENERGY,M(3),X3(3,3),
     &               XDOT3(3,3),CM(10),C11,C12,C19,C20,C24,C25,
     &               NSTEP3,NAME3(3),KZ15,KZ27
      REAL*8  Q1(9),P1(9),Q2(9),P2(9)
*
*
*       Decide path (initialization, termination, KS -> phys, switching).
      IF (KDUM.EQ.4) GO TO 4
      IF (KDUM.GT.1) GO TO 20
*
*       Obtain total kinetic & potential energy.
      ZKE = 0.0D0
      POT = 0.0D0
      DO 3 I = 1,3
          ZKE = ZKE + M(I)*(XDOT3(1,I)**2 + XDOT3(2,I)**2 +
     &                                      XDOT3(3,I)**2)
          DO 2 J = I+1,3
              POT = POT - M(I)*M(J)/SQRT((X3(1,I) - X3(1,J))**2 +
     &                                   (X3(2,I) - X3(2,J))**2 +
     &                                   (X3(3,I) - X3(3,J))**2)
    2     CONTINUE
    3 CONTINUE
*
*       Save twice the initial energy for routine DERQP3.
      ENERGY = ZKE + 2.0D0*POT
      IF (KDUM.EQ.0) GO TO 40
*
*       Introduce physical coordinates & momenta.
    4 DO 6 I = 1,3
          DO 5 K = 1,3
              I1 = 3*I + K - 3
              Q1(I1) = X3(K,I)
              P1(I1) = M(I)*XDOT3(K,I)
    5     CONTINUE
    6 CONTINUE
*
*       Set mass factors for routine DERQP3.
      C11 = 0.25D0/M(1) + 0.25D0/M(3)
      C12 = 0.25D0/M(2) + 0.25D0/M(3)
      C19 = 2.0D0*M(2)*M(3)
      C20 = 2.0D0*M(1)*M(3)
      C24 = 0.25D0/M(3)
      C25 = 2.0D0*M(1)*M(2)
*
*       Define relative coordinates & absolute momenta (eqn (45)).
      DO 7 K = 1,3
          Q2(K) = Q1(K) - Q1(K+6)
          Q2(K+3) = Q1(K+3) - Q1(K+6)
          P2(K+3) = P1(K+3)
    7 CONTINUE
*
*       Expand the variables by relabelling (eqn (46)).
      DO 8 K = 1,3
          Q1(K) = Q2(K)
          Q1(K+4) = Q2(K+3)
          P1(K+4) = P2(K+3)
    8 CONTINUE
*
*       Initialize the redundant variables (eqn (47)).
      Q1(4) = 0.0D0
      Q1(8) = 0.0D0
      P1(4) = 0.0D0
      P1(8) = 0.0D0
*
*       Introduce regularized variables for each KS pair.
      DO 10 KCOMP = 1,2
          K = 4*(KCOMP - 1)
*
*       Form scalar distance from relative separation vector.
          RK = SQRT(Q1(K+1)**2 + Q1(K+2)**2 + Q1(K+3)**2)
*
*       Perform KS coordinate transformation (eqn (48) or (49)).
          IF (Q1(K+1).GE.0.0) THEN
              Q(K+1) = SQRT(0.5D0*(RK + Q1(K+1)))
              Q(K+2) = 0.5D0*Q1(K+2)/Q(K+1)
              Q(K+3) = 0.5D0*Q1(K+3)/Q(K+1)
              Q(K+4) = 0.0D0
          ELSE
              Q(K+2) = SQRT(0.5D0*(RK - Q1(K+1)))
              Q(K+1) = 0.5D0*Q1(K+2)/Q(K+2)
              Q(K+4) = 0.5D0*Q1(K+3)/Q(K+2)
              Q(K+3) = 0.0D0
          END IF
*
*       Set regularized momenta (eqn (50)).
          P(K+1) = 2.0D0*(+Q(K+1)*P1(K+1) + Q(K+2)*P1(K+2) +
     &                                      Q(K+3)*P1(K+3))
          P(K+2) = 2.0D0*(-Q(K+2)*P1(K+1) + Q(K+1)*P1(K+2) +
     &                                      Q(K+4)*P1(K+3))
          P(K+3) = 2.0D0*(-Q(K+3)*P1(K+1) - Q(K+4)*P1(K+2) +
     &                                      Q(K+1)*P1(K+3))
          P(K+4) = 2.0D0*(+Q(K+4)*P1(K+1) - Q(K+3)*P1(K+2) +
     &                                      Q(K+2)*P1(K+3))
   10 CONTINUE
*
      GO TO 40
*
*       Transform each KS pair from regularized to physical variables.
   20 DO 22 KCOMP = 1,2
          K = 4*(KCOMP - 1)
*
*       Obtain relative coordinates (eqn (52)).
          Q1(K+1) = Q(K+1)**2 - Q(K+2)**2 - Q(K+3)**2 + Q(K+4)**2
          Q1(K+2) = 2.0D0*(Q(K+1)*Q(K+2) - Q(K+3)*Q(K+4))
          Q1(K+3) = 2.0D0*(Q(K+1)*Q(K+3) + Q(K+2)*Q(K+4))
*
*       Form product of half Levi-Civita matrix & regularized momenta.
          P1(K+1) = Q(K+1)*P(K+1) - Q(K+2)*P(K+2) - Q(K+3)*P(K+3) +
     &                                              Q(K+4)*P(K+4)
          P1(K+2) = Q(K+2)*P(K+1) + Q(K+1)*P(K+2) - Q(K+4)*P(K+3) -
     &                                              Q(K+3)*P(K+4)
          P1(K+3) = Q(K+3)*P(K+1) + Q(K+4)*P(K+2) + Q(K+1)*P(K+3) +
     &                                              Q(K+2)*P(K+4)
*
*       Evaluate scalar distance.
          RK = Q(K+1)**2 + Q(K+2)**2 + Q(K+3)**2 + Q(K+4)**2
*
*       Set absolute momenta (eqn (53)).
          P1(K+1) = 0.5D0*P1(K+1)/RK
          P1(K+2) = 0.5D0*P1(K+2)/RK
          P1(K+3) = 0.5D0*P1(K+3)/RK
   22 CONTINUE
*
*       Re-define variables in six dimensions (eqn (54)). 
      DO 25 K = 1,3
          Q1(K+3) = Q1(K+4)
          P1(K+3) = P1(K+4)
   25 CONTINUE
*
*       Obtain physical coordinates & momenta in c.m. frame.
      DO 26 K = 1,3
          Q2(K+6) = -(M(1)*Q1(K) + M(2)*Q1(K+3))/(M(1) + M(2) + M(3))
*       Physical coordinates of reference body M(3) (first eqn (55)).
          Q2(K) = Q1(K) + Q2(K+6)
          Q2(K+3) = Q1(K+3) + Q2(K+6)
*       Physical coordinates of body M(1) & M(2) (second eqn (55)).
          P2(K) = P1(K)
          P2(K+3) = P1(K+3)
*       Physical momenta of body M(1) & M(2) (third eqn (55)).
          P2(K+6) = -(P2(K) + P2(K+3))
*       Absolute momentum of reference body M(3) (last eqn (55)).
   26 CONTINUE
*
*       Specify coordinates & velocities in c.m. frame.
      DO 30 I = 1,3
          DO 28 K = 1,3
              I1 = 3*I + K - 3
              X3(K,I) = Q2(I1)
              XDOT3(K,I) = P2(I1)/M(I)
   28     CONTINUE
   30 CONTINUE
*
      IF (KDUM.EQ.3) GO TO 40
*
*       Evaluate total energy & relative energy error.
      S1 = P2(1)**2 + P2(2)**2 + P2(3)**2
      S2 = P2(4)**2 + P2(5)**2 + P2(6)**2
      S3 = P2(7)**2 + P2(8)**2 + P2(9)**2
      ZKE = 0.5D0*(S1/M(1) + S2/M(2) + S3/M(3))
      S1 = M(1)*M(3)/SQRT((Q2(7) - Q2(1))**2 + (Q2(8) - Q2(2))**2 +
     &                                         (Q2(9) - Q2(3))**2)
      S2 = M(2)*M(3)/SQRT((Q2(7) - Q2(4))**2 + (Q2(8) - Q2(5))**2 +
     &                                         (Q2(9) - Q2(6))**2)
      S3 = M(1)*M(2)/SQRT((Q2(4) - Q2(1))**2 + (Q2(5) - Q2(2))**2 +
     &                                         (Q2(6) - Q2(3))**2)
      HT = ZKE - S1 - S2 - S3
*       Current total energy computed from physical variables.
      CM(10) = (HT - 0.5D0*ENERGY)/(0.5D0*ENERGY)
*       Relative energy error with respect to initial value.
*
   40 RETURN
*
      END
