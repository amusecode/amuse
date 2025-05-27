      SUBROUTINE SWEEP(DTCL,RCL)
*
*
*       Enforced KS regularization of wide binaries.
*       --------------------------------------------
*
      INCLUDE 'common6.h'
      SAVE I
*
*
*       Identify wider binaries giving rise to systematic errors.
      ZMX = 100.0*BODYM
      I = IFIRST - 1
    1 I = I + 1
      IF (I.GE.N-1) GO TO 30
      IF (BODY(I).GT.ZMX) GO TO 2
      IF (STEP(I).GT.DTCL.OR.BODY(I).EQ.0.0D0) GO TO 1
    2 CONTINUE
*
*       Search neighbour list of all KS candidates (STEP < DTCL).
      RX2 = 100.0
      NNB1 = LIST(1,I) + 1
      IF (NNB1.EQ.1) GO TO 1
      DO 10 L = 2,NNB1
          J = LIST(L,I)
          RIJ2 = 0.0
          DO 5 K = 1,3
              RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
    5     CONTINUE
          IF (RIJ2.LT.RX2) THEN
              RX2 = RIJ2
              JMIN = J
          END IF
   10 CONTINUE
*       Skip any close c.m. body (small STEP treated by IMPACT).
      IF (JMIN.GT.N) GO TO 1
*
*       Form inverse semi-major axis.
      VIJ2 = 0.0
      DO 15 K = 1,3
          VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,JMIN))**2
   15 CONTINUE
      AINV = 2.0/SQRT(RX2) - VIJ2/(BODY(I) + BODY(JMIN))
*
      IF (AINV.GT.0.0) THEN
          SEMI = 1.0/AINV
          ZMB = BODY(I) + BODY(JMIN)
          EREL = 0.5*ZMB/SEMI
      ELSE
          EREL = 0.0
      END IF
*       Initialize bound KS pairs outside standard parameters.
      IF (AINV.GT.1.0/RCL.OR.EREL.GT.ECLOSE) THEN
          ICOMP = MIN(I,JMIN)
          JCOMP = MAX(I,JMIN)
*       Skip possible case of chain c.m. forming binary.
          IF (NAME(ICOMP).EQ.0.OR.NAME(JCOMP).EQ.0) GO TO 1
*       Ensure most recent velocity used for new KS.
          DO 16 K = 1,3
              X0DOT(K,ICOMP) = XDOT(K,ICOMP)
              X0DOT(K,JCOMP) = XDOT(K,JCOMP)
   16     CONTINUE
          CALL KSREG
          NEWKS = NEWKS + 1
          RI2 = 0.0
          DO 20 K = 1,3
              RI2 = RI2 + (X(K,NTOT) - RDENS(K))**2
   20     CONTINUE
          SEMI = 1.0/AINV
          ECC2 = (1.0 - R(NPAIRS)/SEMI)**2 +
     &                            TDOT2(NPAIRS)**2/(BODY(NTOT)*SEMI)
          ECC = SQRT(ECC2)
          J1 = 2*NPAIRS - 1
          IF (NEWKS.LT.50.OR.EREL.GT.ECLOSE) THEN
             WRITE (6,25)   NAME(J1), NAME(J1+1), LIST(1,J1),
     &                      LIST(1,NTOT), SQRT(RI2), ECC, SEMI,
     &                      GAMMA(NPAIRS)
   25        FORMAT (' ENFORCED KS    NAM NP NNB r E A GAM ',
     &                                2I7,I4,I5,F7.2,F8.3,1P,E10.2,E9.1)
          END IF
      END IF
      IF (I.LT.N-2) GO TO 1
*
   30 RETURN
*
      END
