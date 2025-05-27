      SUBROUTINE NBINT(I,NBFLAG)
*
*
*       Irregular integration.
*       ----------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      REAL*8  XI(3),XIDOT(3),FIRR(3),FREG(3),FD(3),FDUM(3),DV(3)
*
      IF(NBFLAG.EQ.0)THEN
*       Predict current state vector of body #I to order FDOT.
          S = TIME - T0(I)
          DO 11 K = 1,3
           X(K,I) = ((FDOT(K,I)*S + F(K,I))*S + X0DOT(K,I))*S + X0(K,I)
           XDOT(K,I) = (3.0*FDOT(K,I)*S + 2.0*F(K,I))*S + X0DOT(K,I)
   11     CONTINUE
*       Predict coordinates & velocities of neighbours to order FDOT (R.Sp.).
          NNB1 = LIST(1,I) + 1
*
          DO 1 L = 2,NNB1
              J = LIST(L,I)
              S = TIME - T0(J)
              S1 = 1.5*S
              S2 = 2.0*S
              X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
              X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
              X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
              XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
              XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
              XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
   1      CONTINUE
*
      END IF
*         call cputim(tt2)
*         ttpre = ttpre + (tt2-tt1)*60.
*         ttnbp = ttnbp + (tt2-tt1)*60.
*
*       Obtain irregular force & first derivative.
      DO 5 K = 1,3
          XI(K) = X(K,I)
          XIDOT(K) = XDOT(K,I)
          FIRR(K) = 0.0D0
          FD(K) = 0.0D0
    5 CONTINUE
          PHII(I) = 0.D0
          NNB0 = LIST(1,I)
          NIRRF = NIRRF + NNB0
*
*       Assume small mass at centre for special case of no neighbours.
      IF (NNB0.EQ.0) THEN
          RI2 = XI(1)**2 + XI(2)**2 + XI(3)**2
          FIJ = 0.01*BODYM/(RI2*SQRT(RI2))
          RDOT = 3.0*(XI(1)*XIDOT(1) + XI(2)*XIDOT(2) +
     &                                 XI(3)*XIDOT(3))/RI2
          DO 10 K = 1,3
              FIRR(K) = -FIJ*XI(K)
              FD(K) = -(XIDOT(K) - RDOT*XI(K))*FIJ
   10     CONTINUE
          IF (I.GT.N) IPAIR = I - N
          GO TO 70
      END IF
*
*       Choose force loop for single particle or regularized c.m. body.
      IF (I.LE.N) GO TO 20
*
*       Set KS pair index.
      IPAIR = I - N
*       Adopt c.m. approximation for small total perturbation.
      I1 = 2*IPAIR - 1
      IF (LIST(1,I1).GT.0) THEN
*       Obtain irregular force on perturbed c.m. body (including any chain).
          CALL CMFIRR(I,IPAIR,XI,XIDOT,FIRR,FD)
          GO TO 70
      END IF
*
*       Set neighbour number & list index of the last single particle.
   20 NNB1 = NNB0 + 1
      NNB2 = NNB1
   25 IF (LIST(NNB2,I).LE.N) GO TO 30
      NNB2 = NNB2 - 1
      IF (NNB2.GT.1) GO TO 25
*       Include special case of only c.m. neighbours.
      GO TO 40
*
*       Sum over single particles (unperturbed case included).
   30 DO 35 L = 2,NNB2
          K = LIST(L,I)
          A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          DV(1) = XDOT(1,K) - XIDOT(1)
          DV(2) = XDOT(2,K) - XIDOT(2)
          DV(3) = XDOT(3,K) - XIDOT(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
          DR2I = 1.0/RIJ2
          DR3I = BODY(K)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
          FIRR(1) = FIRR(1) + A1*DR3I
          FIRR(2) = FIRR(2) + A2*DR3I
          FIRR(3) = FIRR(3) + A3*DR3I
          PHII(I) = PHII(I) - DR3I*RIJ2
          FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
          FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
          FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
   35 CONTINUE
*
*       See whether any c.m. neighbours should be included.
      IF (NNB2.EQ.NNB1) GO TO 60
*
   40 NNB3 = NNB2 + 1
*       Set index for distinguishing c.m. or resolved components.
      KDUM = 0
*
*       Sum over regularized c.m. neighbours.
      DO 50 L = NNB3,NNB1
          K = LIST(L,I)
          A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          DV(1) = XDOT(1,K) - XIDOT(1)
          DV(2) = XDOT(2,K) - XIDOT(2)
          DV(3) = XDOT(3,K) - XIDOT(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
*       See whether c.m. approximation applies (ignore unperturbed case).
          J = K - N
          KDUM = 2*J - 1
          IF (RIJ2.GT.CMSEP2*R(J)**2.OR.LIST(1,KDUM).EQ.0) GO TO 48
*
          K = KDUM
*       Sum over individual components of pair #J.
   45     A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          DV(1) = XDOT(1,K) - XIDOT(1)
          DV(2) = XDOT(2,K) - XIDOT(2)
          DV(3) = XDOT(3,K) - XIDOT(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
*       Adopt c.m. approximation outside the effective perturber sphere.
   48     DR2I = 1.0/RIJ2
          DR3I = BODY(K)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          FIRR(1) = FIRR(1) + A1*DR3I
          FIRR(2) = FIRR(2) + A2*DR3I
          FIRR(3) = FIRR(3) + A3*DR3I
          FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
          FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
          FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
          IF (K.EQ.KDUM) THEN
              K = K + 1
              GO TO 45
          END IF
   50 CONTINUE
*
*       Include treatment for regularized clump.
   60 IF (NCH.GT.0) THEN
*       Distinguish between chain c.m. and any other particle.
*       Note: in NBODY6++ CHFIRR and FCHAIN are called with IR=1 since
*             chain prediction and perturber list are updated in integrator.
          IF (NAME(I).EQ.0) THEN
              CALL CHFIRR(I,1,XI,XIDOT,FIRR,FD)
          ELSE
*       See if chain perturber list contains body #I.
              NP1 = LISTC(1) + 1
              DO 65 L = 2,NP1
                  J = LISTC(L)
                  IF (J.GT.I) GO TO 70
                  IF (J.EQ.I) THEN
                      CALL FCHAIN(I,1,XI,XIDOT,FIRR,FD)
                      GO TO 70
                  END IF
   65         CONTINUE
          END IF
      END IF
*
*       Check option for external tidal field.
*       Use predicted force of previous step for compatibility with regint.
   70 DT = TIME - T0(I)
      IF (KZ(14).GT.0) THEN
          DO 75 K = 1,3
              FREG(K) = FR(K,I) + DT*FRDOT(K,I)
   75     CONTINUE
          CALL XTRNLF(XI,XIDOT,FIRR,FREG,FD,FDUM,0)
      END IF
*
*       Include the corrector and set new F, FDOT, D1, D2 & D3.
      DTSQ = DT**2
      DT6 = 6.0/(DT*DTSQ)
      DT2 = 2.0/DTSQ
      DTSQ12 = ONE12*DTSQ
      DT13 = ONE3*DT
*
      DO 80 K = 1,3
	  DF = FI(K,I) - FIRR(K)
	  FID = FIDOT(K,I)
	  SUM = FID + FD(K)
	  AT3 = 2.0*DF + DT*SUM
	  BT2 = -3.0*DF - DT*(SUM + FID)
*       Use here new variables for consistency in parallel execution (R.Sp.)
          XN(K,I) = XI(K) + (0.6*AT3 + BT2)*DTSQ12
          XNDOT(K,I) = XIDOT(K) + (0.75*AT3 + BT2)*DT13
*
      FI(K,I) = FIRR(K)
      FIDOT(K,I) = FD(K)
*       Use total force for irregular step (cf. Makino & Aarseth PASJ, 1992).
          FDUM(K) = FIRR(K) + FR(K,I)
*
          D2(K,I) = (3.0*AT3 + BT2)*DT2
          D3(K,I) = AT3*DT6
*       NOTE: These are real derivatives!
   80 continue
*
      TTMP = TSTEP(FDUM,FD,D2(1,I),D3(1,I),ETAI)
      DT0 = TTMP
*
*       Suggestion due to Winston Sweatman
*     DVV = (XDOT(1,I)-X0DOT(1,I))**2 + (XDOT(2,I)-X0DOT(2,I))**2 +
*    &     (XDOT(3,I)-X0DOT(3,I))**2
*     FFD = FIRR(1)**2 + FIRR(2)**2 + FIRR(3)**2
*     ETAIW = ETAI
*     TTMPW = ETAIW*DVV*BODY(I)/FFD
*
*     PRINT*,' irr I=',I,' TTMP,TTMPW,RATIO=',
*    &  TTMP,TTMPW,TTMP/TTMPW
*
*     IF(TTMP.GT.TTMPW)THEN
*     IGT = IGT + 1
*     ELSE
*     ILE = ILE + 1
*     END IF
*     IF(MOD(IGT+ILE,100).EQ.0)PRINT*,' irr IGT,ILE=',IGT,ILE
*
*     TTMP = MAX(TTMPW,TTMP)
*     DT0 = TTMP
*
*     IF (I.GT.N) THEN
*       Check for hierarchical configuration but exclude small perturbations.
*         IF (H(IPAIR).LT.-ECLOSE.AND.KZ(36).GT.0) THEN
*             IF (GAMMA(IPAIR).GT.1.0E-04) THEN
*                 CALL KEPLER(I,TTMP)
*                 DT0 = TTMP
*             END IF
*         END IF
*     END IF
*
*       Include convergence test for large step (cf. Makino, Ap.J. 369, 200).
      IF (TTMP.GT.STEPJ.AND.N.GT.1000) THEN
         DV2 = 0.0
         F2 = 0.0
*       Use only low order predicted value here.
         DO 85 K = 1,3
            DV2 = DV2 + (XIDOT(K) - XNDOT(K,I))**2
            F2 = F2 + FIRR(K)**2
   85    CONTINUE
*       Employ Makino criterion to avoid over-shooting (cf. Book, 2.16).
         DTJ = STEP(I)*(1.0D-06*STEP(I)**2*F2/DV2)**0.1
         TTMP = MIN(TTMP,DTJ)
      END IF
*
*       Select discrete value (increased by 2, decreased by 2 or unchanged).
      IF (TTMP.GT.2.0*STEP(I)) THEN
          IF (DMOD(TIME,2.0*STEP(I)).EQ.0.0D0) THEN
              TTMP = MIN(2.0*STEP(I),1.D0)
          ELSE
              TTMP = STEP(I)
          END IF
      ELSE IF (TTMP.LT.STEP(I)) THEN
          TTMP = 0.5*STEP(I)
            IF (TTMP.GT.DT0) THEN
                TTMP = 0.5*TTMP
          END IF
      ELSE
          TTMP = STEP(I)
      END IF
*
      STEP(I) = TTMP
*
      RETURN
*
      END
