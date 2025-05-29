      SUBROUTINE KSCORR(IPAIR,UI,UIDOT,FP,FD,TD2,TDOT4,TDOT5,TDOT6)
*
*
*       Corrector for KS motion.
*       ------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      REAL*8  UI(4),UIDOT(4),FP(6),FD(6),FREG(4),FRD(4),A1(3,4),A(8),
     &        PR(4),PRD(4),U2(4),U3(4),U4(4),U5(4),FD1(3)
      PARAMETER  (ITLOW=2)
*
*
*       Convert from physical to regularized derivative using T' = R.
      RI = R(IPAIR)
      DO 1 K = 1,3
          FD1(K) = FD(K)
          FD(K) = RI*FD(K)
    1 CONTINUE
*
*       Include KS slow-down factor in the perturbation if ZMOD > 1.
      IF (KZ(26).GT.0) THEN
          IMOD = KSLOW(IPAIR)
          IF (IMOD.GT.1) THEN
              ZMOD = FLOAT(ISLOW(IMOD))
              DO 5 K = 1,3
                  FP(K) = ZMOD*FP(K)
                  FD(K) = ZMOD*FD(K)
                  FD1(K) = ZMOD*FD1(K)
    5         CONTINUE
          END IF
      END IF
*
*       Predict DH from Taylor series derivatives and save temporary value.
      DTU = DTAU(IPAIR)
      DT3 = ONE6*DTU
      DT4 = 0.25D0*DTU
      DH = (((HDOT4(IPAIR)*DT4 + HDOT3(IPAIR))*DT3 +
     &                     0.5D0*HDOT2(IPAIR))*DTU + HDOT(IPAIR))*DTU
*
*       Set time-step factors and copy Stumpff coefficients to scalars.
      DTSQ = DTU**2
      DT6 = 6.0/(DTU*DTSQ)
      DT2 = 2.0/DTSQ
      Z3 = SF(3,IPAIR)
      Z4 = SF(4,IPAIR)
      Z5 = SF(5,IPAIR)
*
*       Define optimized scalars outside iteration loop.
      DT5 = 0.2D0*DTU*Z5
      DTZ = DT4*Z4
      DH0 = (0.5D0*HDOT2(IPAIR)*DTU + HDOT(IPAIR))*DTU
*
*       Use extra iteration and new KSPERT for large perturbation.
*     IF (GAMMA(IPAIR).LT.1.0D-20) THEN
          ITMAX = ITLOW
*     ELSE
*         ITMAX = ITLOW + 1
*         I1 = 2*IPAIR - 1
*         NNB0 = LIST(1,I1)
*         I = N + IPAIR
*         BODYIN = 1.0/BODY(I)
*     END IF
*
*       Perform iteration with or without re-calculating perturbation.
      DO 40 ITER = 1,ITMAX
*
*       Obtain new transformation matrix.
          CALL MATRIX(UI,A1)
*
*       Form twice regularized force and half first derivative of H.
          HD = 0.0D0
          TD2 = 0.0D0
          DO 10 K = 1,4
              A(K) = A1(1,K)*FP(1) + A1(2,K)*FP(2) + A1(3,K)*FP(3)
              A(K+4) = A1(1,K)*FD(1) + A1(2,K)*FD(2) + A1(3,K)*FD(3)
              PR(K) = RI*A(K)
              FREG(K) = DH*UI(K) + PR(K)
              HD = HD + UIDOT(K)*A(K)
              TD2 = TD2 + UI(K)*UIDOT(K)
   10     CONTINUE
*
*       Set regularized velocity matrix (Levi-Civita matrix not required).
          CALL MATRIX(UIDOT,A1)
*
*       Include the whole (L*F)' term in explicit derivatives of FU & H'.
          HD2 = 0.0D0
          DO 15 K = 1,4
              AK4 = A(K+4) + A1(1,K)*FP(1) + A1(2,K)*FP(2) +
     &                                       A1(3,K)*FP(3)
              HD2 = HD2 + (H(IPAIR)*UI(K) + FREG(K))*A(K) +
     &                                      2.0D0*UIDOT(K)*AK4
              PRD(K) = 0.5D0*RI*AK4 + TD2*A(K)
              FRD(K) = HD*UI(K) + 0.5D0*DH*UIDOT(K) + PRD(K)
*       Form the regularized perturbation and modified force.
              PR(K) = 0.5D0*PR(K)
              FREG(K) = 0.5D0*FREG(K)
   15     CONTINUE
*
*       Determine new derivatives of U evaluated at beginning of the step.
          DO 20 K = 1,4
              DF = FP0(K,IPAIR) - FREG(K)
              SUM = FD0(K,IPAIR) + FRD(K)
	      BT2 = -3.0D0*DF - (SUM + FD0(K,IPAIR))*DTU
	      AT3 = 2.0D0*DF + SUM*DTU
*
              U2(K) = 0.5D0*H0(IPAIR)*U0(K,IPAIR) + FP0(K,IPAIR)
              U3(K) = 0.5D0*H0(IPAIR)*UDOT(K,IPAIR) + FD0(K,IPAIR)
              U4(K) = 0.5D0*H0(IPAIR)*U2(K) + BT2*DT2
              U5(K) = 0.5D0*H0(IPAIR)*U3(K) + AT3*DT6
   20     CONTINUE
*
*       Improve the solution of U and UDOT.
          DO 25 K = 1,4
              UI(K) = ((((U5(K)*DT5 + U4(K)*Z4)*DT4 + U3(K))*DT3 +
     &                 0.5*U2(K))*DTU + UDOT(K,IPAIR))*DTU + U0(K,IPAIR)
              UIDOT(K) = (((U5(K)*DTZ + U4(K)*Z3)*DT3 +
     &                       0.5*U3(K))*DTU + U2(K))*DTU + UDOT(K,IPAIR)
   25     CONTINUE
*
*       Update the physical distance for next iteration or final solution.
          RI = UI(1)**2 + UI(2)**2 + UI(3)**2 + UI(4)**2
*
*       Choose between old and new perturbation (skip last time).
*         IF (ITER.GT.ITMAX) THEN
*         IF (GAMMA(IPAIR).LT.1.0D-04) THEN
*             DO 30 K = 1,3
*                 FD(K) = RI*FD1(K)
*  30         CONTINUE
*         ELSE
*       Transform to improved coordinates & velocities.
*             CALL KSTRAN(IPAIR,I1,I,BODYIN,RI,UI,UIDOT,XI,VI)
*
*       Re-calculate the perturbing force & derivative.
*             CALL KSPERT(I1,NNB0,XI,VI,FP,FD)
*
*       Convert to regularized derivative (note slow-down not active).
*             DO 35  K = 1,3
*                 FD(K) = RI*FD(K)
*  35         CONTINUE
*         END IF
*         END IF
*
*       Re-evaluate DH by adding Hermite corrector.
          HD = 2.0D0*HD
          DHD = HDOT(IPAIR) - HD
          SUM = HDOT2(IPAIR) + HD2
          BT2 = -3.0D0*DHD - (SUM + HDOT2(IPAIR))*DTU
          AT3 = 2.0D0*DHD + SUM*DTU
          DH = DH0 + (0.25D0*AT3 + ONE3*BT2)*DTU
   40 CONTINUE
*
*       Copy final values and set higher derivatives.
      H(IPAIR) = H(IPAIR) + DH
      DO 50 K = 1,4
          U(K,IPAIR) = UI(K)
          U0(K,IPAIR) = UI(K)
          UDOT(K,IPAIR) = UIDOT(K)
          FUDOT2(K,IPAIR) = U4(K) + U5(K)*DTU
          FUDOT3(K,IPAIR) = U5(K)
          AK4 = A(K+4) + A1(1,K)*FP(1) + A1(2,K)*FP(2) + A1(3,K)*FP(3)
          PR(K) = 0.5D0*RI*A(K)
          PRD(K) = 0.5D0*RI*AK4 + TD2*A(K) + 0.5D0*HD*UI(K)
          FP0(K,IPAIR) = PR(K)
          FD0(K,IPAIR) = PRD(K)
   50 CONTINUE
*
*       Save new derivatives of H.
      HDOT(IPAIR) = HD
      HDOT2(IPAIR) = HD2
      HDOT3(IPAIR) = (3.0D0*AT3 + BT2)*DT2
      HDOT4(IPAIR) = AT3*DT6
*
*       Set new FU/2 & FUDOT/6 and form scalar terms for time derivatives.
      TD2 = 0.0D0
      TD3 = 0.0D0
      TDOT4 = 0.0D0
      TDOT5 = 0.0D0
      TDOT6 = 0.0D0
      H2 = 0.5D0*H(IPAIR)
      DO 60 K = 1,4
          U2K = PR(K) + H2*UI(K)
          U3K = PRD(K) + H2*UIDOT(K)
          FU(K,IPAIR) = 0.5D0*U2K
          FUDOT(K,IPAIR) = ONE6*U3K
          TD2 = TD2 + UI(K)*UIDOT(K)
          TD3 = TD3 + UIDOT(K)**2 + UI(K)*U2K
          TDOT4 = TDOT4 + UI(K)*U3K + 3.0D0*UIDOT(K)*U2K
          TDOT5 = TDOT5 + 0.5D0*FUDOT2(K,IPAIR)*UI(K) +
     &                    2.0D0*U3K*UIDOT(K) + 1.5D0*U2K**2
          TDOT6 = TDOT6 + U5(K)*UI(K) + 5.0D0*U4(K)*UIDOT(K) +
     &                    10.0D0*U3K*U2K
*       Note that TDOT4/2, TDOT5/4 and TDOT6/2 are calculated.
   60 CONTINUE
*
*       Save distance and second & third time derivatives.
      R(IPAIR) = RI
      TDOT2(IPAIR) = 2.0D0*TD2
      TDOT3(IPAIR) = 2.0D0*TD3
*
      RETURN
*
      END
