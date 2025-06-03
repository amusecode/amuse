      SUBROUTINE RESOLV(IPAIR,KDUM)
*
*
*       Transformation of KS variables.
*       -------------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      REAL*8  Q(3),RDOT(3),UI(4),V(4),A1(3,4)
*
*
*       Note: KDUM = 1 for prediction of U & UDOT and = 2 for termination.
      I1 = 2*IPAIR - 1
      I2 = 2*IPAIR
*
*       Ensure appropriate prediction (second pair in binary collision).
      IF (LIST(1,I1).EQ.0.OR.T0(I1).EQ.TIME) THEN
*       Set current values of regularized coordinates & velocities.
          DO 1 K = 1,4
              UI(K) = U0(K,IPAIR)
              V(K) = UDOT(K,IPAIR)
    1     CONTINUE
*
*       Copy binding energy for routine ENERGY.
          HT = H(IPAIR)
          GO TO 4
      END IF
*
*       Predict U & UDOT to order FUDOT3 using third-order time inversion.
      A2 = 1.0/R(IPAIR)
      A3 = A2*(TIME - T0(I1))
*
*       See whether the time interval should be modified by KSLOW procedure.
      IF (KSLOW(IPAIR).GT.1) THEN
          IMOD = KSLOW(IPAIR)
          A3 = A3/FLOAT(ISLOW(IMOD))
      END IF
*
*       Expand regularized interval to third order.
      A4 = 3.0D0*TDOT2(IPAIR)**2*A2 - TDOT3(IPAIR)
      DTU = ((ONE6*A4*A3 - 0.5D0*TDOT2(IPAIR))*A2*A3 + 1.0)*A3
*       Apply safety test near small pericentre or for unperturbed motion.
      IF (DTU.GT.DTAU(IPAIR)) DTU = 0.8*DTAU(IPAIR)
      IF (DTU.LT.0.0) DTU = 0.0
*
      DTU1 = 0.2D0*DTU
      DTU2 = DTU/24.0D0
      DO 3 K = 1,4
          UI(K) = ((((FUDOT3(K,IPAIR)*DTU1 + FUDOT2(K,IPAIR))*DTU2 +
     &                          FUDOT(K,IPAIR))*DTU + FU(K,IPAIR))*DTU +
     &                                  UDOT(K,IPAIR))*DTU + U0(K,IPAIR)
          V(K) = (((FUDOT3(K,IPAIR)*DTU2 + ONE6*FUDOT2(K,IPAIR))*DTU +
     &              3.0D0*FUDOT(K,IPAIR))*DTU + 2.0D0*FU(K,IPAIR))*DTU +
     &                                                     UDOT(K,IPAIR)
    3 CONTINUE
*
*       Predict current binding energy per unit mass for routine ADJUST.
      HT = (((HDOT4(IPAIR)*DTU2 + ONE6*HDOT3(IPAIR))*DTU +
     &             0.5D0*HDOT2(IPAIR))*DTU + HDOT(IPAIR))*DTU + H(IPAIR)
*
*       Form current transformation matrix and two-body separation.
    4 CALL MATRIX(UI,A1)
      RI = UI(1)**2 + UI(2)**2 + UI(3)**2 + UI(4)**2
*
*       Obtain relative coordinates & velocities.
      DO 7 J = 1,3
          Q(J) = 0.0D0
          RDOT(J) = 0.0D0
          DO 6 K = 1,4
              Q(J) = Q(J) + A1(J,K)*UI(K)
              RDOT(J) = RDOT(J) + 2.0D0*A1(J,K)*V(K)/RI
    6     CONTINUE
    7 CONTINUE
*
      I = N + IPAIR
      IF (KDUM.EQ.1) GO TO 12
*
*       Predict improved X & X0DOT for c.m. & components at termination.
      DT = TIME - T0(I)
      A0 = 0.05*DT
      A2 = 0.25*DT
      A3 = T0(I) - T0R(I)
      A4 = 0.5*A3
      T0(I) = TIME
*
      DO 10 K = 1,3
          FK = ((ONE6*D3R(K,I)*A3 + 0.5*D2R(K,I))*A3 + D1R(K,I))*A3
     &                                              + D0R(K,I) + D0(K,I)
          F1DOTK = (D3R(K,I)*A4 + D2R(K,I))*A3 + D1R(K,I) + D1(K,I)
          F2DOTK = 0.5*(D3R(K,I)*A3 + D2R(K,I) + D2(K,I))
          F3DOTK = ONE6*(D3R(K,I) + D3(K,I))
          X(K,I) = ((((F3DOTK*A0 + ONE12*F2DOTK)*DT +
     &                   ONE6*F1DOTK)*DT + 0.5*FK)*DT + X0DOT(K,I))*DT +
     &                   X0(K,I)
*       Note: Update X0 for call from KSTERM after RESOLV call in MDOT.
          X0(K,I) = X(K,I)
          X0DOT(K,I)  = (((F3DOTK*A2 + ONE3*F2DOTK)*DT +
     &                            0.5*F1DOTK)*DT + FK)*DT + X0DOT(K,I)
          XDOT(K,I) = X0DOT(K,I)
          X0DOT(K,I1) = X0DOT(K,I) + BODY(I2)*RDOT(K)/BODY(I)
          X0DOT(K,I2) = X0DOT(K,I) - BODY(I1)*RDOT(K)/BODY(I)
   10 CONTINUE
*
*       Set coordinates & velocities (XDOT for KDUM = 1) of KS components.
   12 DO 15 K = 1,3
          X(K,I1) = X(K,I) + BODY(I2)*Q(K)/BODY(I)
          X(K,I2) = X(K,I) - BODY(I1)*Q(K)/BODY(I)
          XDOT(K,I1) = XDOT(K,I) + BODY(I2)*RDOT(K)/BODY(I)
          XDOT(K,I2) = XDOT(K,I) - BODY(I1)*RDOT(K)/BODY(I)
   15 CONTINUE
*
*       Check updating of KS variables (KDUM = 3 in routine MDOT).
      IF (KDUM.LT.3) GO TO 20
*
      DO 18 K = 1,4
          U(K,IPAIR) = UI(K)
          U0(K,IPAIR) = UI(K)
          UDOT(K,IPAIR) = V(K)
   18 CONTINUE
      R(IPAIR) = RI
      H(IPAIR) = HT
      T0(I1) = TIME
      T0(I) = TIME
*
   20 RETURN
*
      END
