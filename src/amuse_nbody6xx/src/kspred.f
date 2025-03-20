      SUBROUTINE KSPRED(IPAIR,I1,I,BODYIN,UI,UIDOT,XI,VI)
*
*
*       Prediction for KS regularization.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (ONE24=1.0/24.0D0,ONE120=1.0/120.0D0)
      REAL*8  UI(4),UIDOT(4),XI(6),VI(6),RDOT(3),A1(3,4)
*
*
*       Add body #I to the perturber list for prediction.
      NNB2 = LIST(1,I1) + 2
      LIST(NNB2,I1) = I
*
*       Predict coordinates & velocities of perturbers & c.m. to order FDOT.
      DO 10 L = 2,NNB2
          J = LIST(L,I1)
          S = TIME - T0(J)
          S1 = 1.5*S
          S2 = 2.0*S
          X(1,J) = ((FDOT(1,J)*S + F(1,J))*S + X0DOT(1,J))*S + X0(1,J)
          X(2,J) = ((FDOT(2,J)*S + F(2,J))*S + X0DOT(2,J))*S + X0(2,J)
          X(3,J) = ((FDOT(3,J)*S + F(3,J))*S + X0DOT(3,J))*S + X0(3,J)
	  XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
	  XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
	  XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
   10 CONTINUE
*
*       Copy Stumpff coefficients to scalars.
      Z3 = SF(3,IPAIR)
      Z4 = SF(4,IPAIR)
      Z5 = SF(5,IPAIR)
*
*       Predict KS motion to full accuracy including Stumpff functions.
      DTU = DTAU(IPAIR)
      DO 20 K = 1,4
          U2 = FU(K,IPAIR)
          U3 = FUDOT(K,IPAIR)
          U4 = ONE24*FUDOT2(K,IPAIR)
          U5 = ONE120*FUDOT3(K,IPAIR)
*
          UI(K) = ((((U5*Z5*DTU + U4*Z4)*DTU + U3)*DTU + U2)*DTU +
     &                                 UDOT(K,IPAIR))*DTU + U0(K,IPAIR)
          UIDOT(K) = (((5.0*U5*Z4*DTU + 4.0D0*U4*Z3)*DTU + 3.0D0*U3)*DTU
     &                                     + 2.0*U2)*DTU + UDOT(K,IPAIR)
   20 CONTINUE
*
      R(IPAIR) = UI(1)**2 + UI(2)**2 + UI(3)**2 + UI(4)**2
*
*       Form relative coordinates obtained from explicit KS transformation.
      Q1 = UI(1)**2 - UI(2)**2 - UI(3)**2 + UI(4)**2
      Q2 = UI(1)*UI(2) - UI(3)*UI(4)
      Q3 = UI(1)*UI(3) + UI(2)*UI(4)
      Q2 = Q2 + Q2
      Q3 = Q3 + Q3
*
*       Assign global coordinates of regularized components.
      A2 = BODY(I1+1)*BODYIN
      XI(1) = X(1,I) + A2*Q1
      XI(2) = X(2,I) + A2*Q2
      XI(3) = X(3,I) + A2*Q3
      XI(4) = XI(1) - Q1
      XI(5) = XI(2) - Q2
      XI(6) = XI(3) - Q3
*
*       Set current transformation matrix.
      CALL MATRIX(UI,A1)
*
*       Obtain relative velocities from KS transformation.
      RINV = 2.0D0/R(IPAIR)
      DO 30 L = 1,3
          RDOT(L) = 0.0D0
          DO 25 K = 1,4
              RDOT(L) = RDOT(L) + A1(L,K)*UIDOT(K)
   25     CONTINUE
          RDOT(L) = RDOT(L)*RINV
   30 CONTINUE
*
*       Set global velocities of KS components.
      VI(1) = XDOT(1,I) + A2*RDOT(1)
      VI(2) = XDOT(2,I) + A2*RDOT(2)
      VI(3) = XDOT(3,I) + A2*RDOT(3)
      VI(4) = VI(1) - RDOT(1)
      VI(5) = VI(2) - RDOT(2)
      VI(6) = VI(3) - RDOT(3)
*
      RETURN
*
      END
