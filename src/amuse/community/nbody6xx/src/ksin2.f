      SUBROUTINE KSIN2(ICASE)
*
*
*       Initialization of hierarchical KS.
*       ----------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  Q(3),RDOT(3),UI(4),VI(4),A1(3,4)
      CHARACTER*8  WHICH1
*
*
*       Copy KS pair index and set c.m. index.
      IPAIR = KSPAIR
      ICM = N + IPAIR
*
*       Distinguish first or second stage (termination: ICASE > 1 only).
      IF (ICASE.GT.1) GO TO 50
*
*       Specify mass & name for new c.m. and initialize radius & type.
      BODYI = BODY(ICOMP)
      BODY(ICM) = BODY(ICOMP) + BODY(JCOMP)
      NAME(ICM) = NZERO + NAME(ICOMP)
      RADIUS(ICM) = 0.0
      TEV(ICM) = 1.0E+10
      TEV0(ICM) = 1.0E+10
      BODY0(ICM) = BODYM
      EPOCH(ICM) = TIME*TSTAR
      KSTAR(ICM) = 0
*
*       Define relative coordinates and velocities in physical units.
      DO 20 K = 1,3
          Q(K) = X(K,ICOMP) - X(K,JCOMP)
          RDOT(K) = X0DOT(K,ICOMP) - X0DOT(K,JCOMP)
   20 CONTINUE
*
*       Introduce regularized variables using definition of 1985 paper.
      R(IPAIR) = SQRT(Q(1)**2 + Q(2)**2 + Q(3)**2)
*
*       Initialize the regularized coordinates according to sign of Q(1).
      IF (Q(1).LE.0.0D0) THEN
          UI(3) = 0.0D0
          UI(2) = SQRT(0.5D0*(R(IPAIR) - Q(1)))
          UI(1) = 0.5D0*Q(2)/UI(2)
          UI(4) = 0.5D0*Q(3)/UI(2)
      ELSE
          UI(4) = 0.0D0
          UI(1) = SQRT(0.5D0*(R(IPAIR) + Q(1)))
          UI(2) = 0.5D0*Q(2)/UI(1)
          UI(3) = 0.5D0*Q(3)/UI(1)
      END IF
*
*       Set current transformation matrix.
      CALL MATRIX(UI,A1)
*
*       Form regularized velocity and set initial KS coordinates & TDOT2.
      TDOT2(IPAIR) = 0.0D0
      DO 30 K = 1,4
          UDOT(K,IPAIR) = 0.50D0*(A1(1,K)*RDOT(1) + A1(2,K)*RDOT(2) +
     &                                                  A1(3,K)*RDOT(3))
*       Note that A1(J,K) is the transpose of A1(K,J).
          U(K,IPAIR) = UI(K)
          U0(K,IPAIR) = U(K,IPAIR)
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0D0*UI(K)*UDOT(K,IPAIR)
   30 CONTINUE
*
*       Evaluate initial binding energy per unit mass and initialize T0.
      H(IPAIR) = (2.0D0*(UDOT(1,IPAIR)**2 + UDOT(2,IPAIR)**2 +
     &                   UDOT(3,IPAIR)**2 + UDOT(4,IPAIR)**2) -
     &                                              BODY(ICM))/R(IPAIR)
      T0(2*IPAIR-1) = TIME
*
*       Copy old c.m. mass to first KS component.
      BODY(2*IPAIR-1) = BODYI
*
*       Define c.m. coordinates & velocities.
      DO 40 K = 1,3
          X(K,ICM) = (BODYI*X(K,ICOMP) + BODY(JCOMP)*X(K,JCOMP))/
     &                                                        BODY(ICM)
          X0DOT(K,ICM) = (BODYI*X0DOT(K,ICOMP) + BODY(JCOMP)*
     &                                        X0DOT(K,JCOMP))/BODY(ICM)
          XDOT(K,ICM) = X0DOT(K,ICM)
          X0(K,ICM) = X(K,ICM)
   40 CONTINUE
*
*       Exit after initialization.
      GO TO 100
*
*       Check for sufficient perturbers (RMAX = apocentre set in IMPACT).
   50 NNB = LIST(1,ICM)
      IF (SQRT(CMSEP2)*RMAX.GT.2.0*RS(ICM)) THEN
          FAC = FLOAT(NNBMAX)/FLOAT(NNB)
          RSI = FAC**0.3333*RS(ICM)
          IF (RSI.GT.RS(ICM)) THEN
              CALL NBLIST(ICM,RSI)
          END IF
      END IF
*
*       Initialize force polynomial for c.m. using new neighbour list.
      CALL FPOLY1(ICM,ICM,0)
      CALL FPOLY2(ICM,ICM,0)
*
*       Form perturber list.
      CALL KSLIST(IPAIR)
*
*       Transform any unperturbed hard binary to apocentre and set time-step.
      IMOD = 1
      EB = H(IPAIR)*BODY(ICOMP)*BODY(JCOMP)/BODY(ICM)
*       Suppress the following for now.
      IF (LIST(1,ICOMP).LT.0.AND.EB.LT.EBH) THEN
          SEMI = -0.5*BODY(ICM)/H(IPAIR)
          TK = TWOPI*ABS(SEMI)*SQRT(ABS(SEMI)/BODY(ICM))
          IF (IPHASE.NE.7) THEN
              DO 55 K = 1,4
                  UI(K) = U(K,IPAIR)
                  VI(K) = UDOT(K,IPAIR)
   55         CONTINUE
*       Determine pericentre time (TP < 0 if TDOT2 < 0) and add TK/2.
              CALL TPERI(SEMI,UI,VI,BODY(ICM),TP)
              STEP(ICOMP) = 0.5*MIN(TK,STEP(ICM)) - TP
*       Transform KS variables to peri and by pi/2 to apocentre (skip apo).
              IF (ABS(TDOT2(IPAIR)).GT.1.0E-12.OR.R(IPAIR).LT.SEMI) THEN
                  CALL KSPERI(IPAIR)
                  CALL KSAPO(IPAIR)
              ELSE IF (TDOT2(IPAIR).GT.0.0) THEN
                  TDOT2(IPAIR) = -1.0E-20
              END IF
          END IF
*       Estimate an appropriate KS slow-down index for G < GMIN.
          IF (KZ(26).GT.0.AND.STEP(ICM).GT.TK) THEN
              IMOD = 1 + LOG(STEP(ICM)/TK)/0.69
              IMOD = MIN(IMOD,5)
          END IF
      END IF
*
*       Specify zero membership and large step for second component.
      LIST(1,JCOMP) = 0
      STEP(JCOMP) = 1.0E+06
      STEPR(JCOMP) = 1.0E+06
*
*       Obtain polynomials for perturbed KS motion (standard case).
      CALL KSPOLY(IPAIR,IMOD)
*
*       Increase regularization counters (NKSHYP for hyperbolic orbits).
      NKSREG = NKSREG + 1
*
      IF (KZ(10).GT.0.AND.ICASE.GE.2) THEN
          RI = SQRT((X(1,ICM) - RDENS(1))**2 +
     &              (X(2,ICM) - RDENS(2))**2 +
     &              (X(3,ICM) - RDENS(3))**2)
          WHICH1 = ' MERGE2 '
          IF (ICASE.EQ.3) WHICH1 = ' RESET2 '
          if(rank.eq.0)
     &    WRITE (6,60)  WHICH1, TIME+TOFF, NAME(ICOMP), NAME(JCOMP),
     &                  DTAU(IPAIR), R(IPAIR), RI, H(IPAIR), IPAIR,
     &                  GAMMA(IPAIR), STEP(ICM), LIST(1,ICOMP)
   60     FORMAT (/,' NEW',A8,'   T =',F8.2,2I6,F12.3,1P,E10.1,0P,
     &                            F8.2,F9.2,I5,F8.3,1P,E10.1,0P,I5)
      END IF
*
*       See whether either component has been regularized recently.
      NNB = LISTD(1) + 1
      K = 0
*       Check case of initial binary and loop over disrupted pairs.
      IF (IABS(NAME(ICOMP) - NAME(JCOMP)).EQ.1) K = -1
      DO 80 L = 2,NNB
          IF (NAME(ICOMP).EQ.LISTD(L).OR.NAME(JCOMP).EQ.LISTD(L)) K = -1
   80 CONTINUE
*
*       Ensure that mergers are treated as new binaries.
      IF (IPHASE.EQ.6) K = 0
*       Set flags to distinguish primordial binaries & standard KS motion.
      LIST(2,JCOMP) = K
      KSLOW(IPAIR) = 1
*
*       Check diagnostic output of new hard binary.
      IF (KZ(8).GT.0.AND.K.EQ.0) THEN
          IF (EB.GT.EBH) GO TO 100
          SEMI = -0.5*BODY(ICM)/H(IPAIR)
          RI = SQRT((X(1,ICM) - RDENS(1))**2 +
     &              (X(2,ICM) - RDENS(2))**2 +
     &              (X(3,ICM) - RDENS(3))**2)
          if(rank.eq.0)
     &    WRITE (8,90)  TIME+TOFF, NAME(ICOMP), NAME(JCOMP), K,
     &                  BODY(ICOMP), BODY(JCOMP), EB, SEMI, R(IPAIR),
     &                  GAMMA(IPAIR), RI
   90     FORMAT (' NEW BINARY   T =',F7.1,'  NAME = ',2I6,I3,
     &                        '  M =',2F8.4,'  EB =',F9.4,'  A =',F8.5,
     &                          '  R =',F8.5,'  G =',F6.3,'  RI =',F5.2)
      END IF
*
  100 RETURN
*
      END
