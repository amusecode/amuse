      SUBROUTINE KSRES2(J,J1,J2,RIJ2)
*
*
*       Coordinates & velocities of KS pair.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      REAL*8  UI(4),RDOT(3),V(4),A1(3,4)
*
*
*       Skip merged ghost binary (defined by zero total mass).
      IF (BODY(N+J).LE.0.0D0) GO TO 40
*
*       Resolve components of pair #J at curent time.
      J2 = J + J
      J1 = J2 - 1
      A2 = 1.0/R(J)
      A3 = A2*(TIME - T0(J1))
      IF (KSLOW(J).GT.1) THEN
          IMOD = KSLOW(J)
          A3 = A3/FLOAT(ISLOW(IMOD))
      END IF
*
*       Decide appropriate order for interpolation & prediction.
      IF (RIJ2.GT.625.0*R(J)**2) THEN
*       Convert physical interval to regularized time (second order only).
          DTU = (1.0 - 0.5D0*TDOT2(J)*A2*A3)*A3
          IF (ABS(DTU).GT.DTAU(J)) DTU = DTAU(J)
*
*       Predict regularized coordinates of distant pair to second order.
          DO 10 K = 1,4
              UI(K) = (FU(K,J)*DTU + UDOT(K,J))*DTU + U0(K,J)
              V(K) = (3.0*FUDOT(K,J)*DTU + 2.0*FU(K,J))*DTU + UDOT(K,J)
   10     CONTINUE
      ELSE
*        Expand regularized time interval to third order.
          A4 = 3.0D0*TDOT2(J)**2*A2 - TDOT3(J)
          DTU = ((ONE6*A4*A3 - 0.5D0*TDOT2(J))*A2*A3 + 1.0)*A3
*       Apply safety test near small pericentre or for unperturbed motion.
          IF (DTU.GT.DTAU(J)) DTU = 0.8*DTAU(J)
          IF (LIST(1,J1).EQ.0) DTU = 0.0
*
*       Predict regularized coordinates to third order.
          DTU1 = DTU/24.0D0
          DTU2 = ONE6*DTU
          DO 20 K = 1,4
              UI(K) = (((FUDOT2(K,J)*DTU1 + FUDOT(K,J))*DTU +
     &                           FU(K,J))*DTU + UDOT(K,J))*DTU + U0(K,J)
              V(K) = ((FUDOT2(K,J)*DTU2 + 3.0D0*FUDOT(K,J))*DTU +
     &                                    2.0D0*FU(K,J))*DTU + UDOT(K,J)
   20     CONTINUE
      END IF
*
*       Employ KS transformation.
      Q1 = UI(1)**2 - UI(2)**2 - UI(3)**2 + UI(4)**2
      Q2 = UI(1)*UI(2) - UI(3)*UI(4)
      Q3 = UI(1)*UI(3) + UI(2)*UI(4)
      Q2 = Q2 + Q2
      Q3 = Q3 + Q3
      J3 = N + J
      A2 = BODY(J2)/BODY(J3)
*
*       Set global coordinates of regularized components.
      X(1,J1) = X(1,J3) + A2*Q1
      X(2,J1) = X(2,J3) + A2*Q2
      X(3,J1) = X(3,J3) + A2*Q3
      X(1,J2) = X(1,J1) - Q1
      X(2,J2) = X(2,J1) - Q2
      X(3,J2) = X(3,J1) - Q3
*
*       Set current transformation matrix and two-body separation.
      CALL MATRIX(UI,A1)
      RI = UI(1)**2 + UI(2)**2 + UI(3)**2 + UI(4)**2
      RINV = 2.0/RI
*
*       Obtain relative velocities from KS transformation.
      DO 30 L = 1,3
          RDOT(L) = 0.0D0
          DO 25 K = 1,4
              RDOT(L) = RDOT(L) + A1(L,K)*V(K)*RINV
   25     CONTINUE
   30 CONTINUE
*
*       Set global velocities of KS components.
      XDOT(1,J1) = XDOT(1,J3) + A2*RDOT(1)
      XDOT(2,J1) = XDOT(2,J3) + A2*RDOT(2)
      XDOT(3,J1) = XDOT(3,J3) + A2*RDOT(3)
      XDOT(1,J2) = XDOT(1,J1) - RDOT(1)
      XDOT(2,J2) = XDOT(2,J1) - RDOT(2)
      XDOT(3,J2) = XDOT(3,J1) - RDOT(3)
*
   40 RETURN
*
      END
