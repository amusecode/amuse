      SUBROUTINE XCPRED(KCASE)
*
*
*       Prediction of global chain coordinates.
*       ---------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  M,MASS,MC,MIJ,MKK
      PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
*
*
*       Check indicator for prediction of perturbers & c.m.
      IF (KCASE.EQ.0) GO TO 4
*
*       Check adding chain c.m. #ICH to perturber list for prediction.
      IF (KCASE.EQ.1) THEN
          NNB2 = LISTC(1) + 2
          LISTC(NNB2) = ICH
      ELSE
*      Note KCASE = 2 for chain c.m. prediction at new block-time.
          NNB2 = LISTC(1) + 1
      END IF
*
*       Predict coordinates of perturbers & possibly c.m. to order FDOT.
      DO 1 L = 2,NNB2
          J = LISTC(L)
          S = TIME - T0(J)
*       Do not allow prediction outside range (NB! No bad effects in DIFSY1).
*         S = MIN(S,STEP(J))
          X(1,J) = ((FDOT(1,J)*S + F(1,J))*S + X0DOT(1,J))*S + X0(1,J)
          X(2,J) = ((FDOT(2,J)*S + F(2,J))*S + X0DOT(2,J))*S + X0(2,J)
          X(3,J) = ((FDOT(3,J)*S + F(3,J))*S + X0DOT(3,J))*S + X0(3,J)
*       Note most recent velocities are used for perturber integration.
    1 CONTINUE
*
*       Obtain global coordinates & velocities from current chain & c.m.
   4  LK = 0
      DO 10 I = 1,NN
          DO 5 K = 1,3
              LK = LK + 1
              XC(K,I) = XCH(LK) + X(K,ICH)
              UC(K,I) = VCH(LK) + XDOT(K,ICH)
    5     CONTINUE
   10 CONTINUE
*
      RETURN
*
      END
