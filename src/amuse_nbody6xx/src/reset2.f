      SUBROUTINE RESET2
*
*
*      Termination of double hierarchy. 
*      --------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      CHARACTER*11  WHICH1
      REAL*8  XX(3,3),VV(3,3)
      SAVE ISTAB
      DATA ISTAB /0/
*
*
*       Set index of disrupted pair and save output parameters.
      IPAIR = KSPAIR
      I = N + IPAIR
      E1 = BODY(2*IPAIR-1)*BODY(2*IPAIR)*H(IPAIR)/BODY(I)
      G1 = GAMMA(IPAIR)
      R1 = R(IPAIR)
      SEMI1 = -0.5*BODY(I)/H(IPAIR)
      ECC2 = (1.0 - R1/SEMI1)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI1)
      ECC = SQRT(ECC2)
*
*       Locate current position in the merger table.
      IMERGE = 0
      DO 1 K = 1,NMERGE
          IF (NAMEM(K).EQ.NAME(I)) IMERGE = K
    1 CONTINUE
*
*       Produce info on quintuplet or sextuplet (NAMEG has outer c.m. name).
      IF (NAMEG(IMERGE).GT.NZERO) THEN
          RI2 = 0.0
          DO 2 K = 1,3
              RI2 = RI2 + (X(K,I) - RDENS(K))**2
    2     CONTINUE
          RI = SQRT(RI2)
          PMIN = SEMI1*(1.0 - ECC)
          EB = 0.5*BODY(2*IPAIR-1)*BODY(2*IPAIR)/SEMI1
          EB = EB*FLOAT(N-NPAIRS)/ZKIN
          ZM = BODY(I)*SMU
          WHICH1 = ' QUINTUPLET'
*       Determine previous merger index with NAME greater by 2*NZERO.
          JM = 1
          DO 3 K = 1,NMERGE
              IF (NAMEM(K).EQ.NAME(I) + 2*NZERO) JM = K
    3     CONTINUE
*       Use previous merger index to identify ghost from earlier level.
          DO 4 J = IFIRST,NTOT
              IF (NAME(J).EQ.NAMEG(JM)) JG = J
    4     CONTINUE
          IF (JG.GT.N) WHICH1 = ' SEXTUPLET '
          if(rank.eq.0)
     &    WRITE (6,5)  WHICH1, TIME+TOFF, ZM, NAME(2*IPAIR-1),
     &                 NAME(2*IPAIR), NAMEG(IMERGE), RI, ECC, EB, SEMI1,
     &                 PMIN, G1
    5     FORMAT (' END',A11,'   T MT NM1 NM2 NM3 RI E1 EB1 A1 PM G1 ',
     &                           F9.2,F6.2,3I6,2F6.2,F6.1,1P,3E10.2)
      END IF
*
*       Check presence of [[B,S],S] quartet or [[B,B],S] quintuplet.
      IF (NAMEG(IMERGE).GT.0.AND.NAMEG(IMERGE).LE.NZERO) THEN
          I1 = 2*IPAIR - 1
          I2 = 2*IPAIR
          CALL FINDJ(I1,JG,IM)
          ZM = BODY(I)*SMU
          IF (NAME(I2).LE.NZERO) THEN
              if(rank.eq.0)
     &        WRITE (6,9)  TIME+TOFF, ZM, NAME(I1), NAME(I2), NAME(JG),
     &                     ECC, SEMI1, R1, G1
    9         FORMAT (/,' END QUARTET    T MT NM1 NMG NM3 E1 A1 R1 G1 ',
     &                                   F9.2,F6.2,3I6,F6.2,1P,3E10.2)
          ELSE
              if(rank.eq.0)
     &        WRITE (6,11)  TIME+TOFF, ZM, NAME(I1), NAME(I2), NAME(JG),
     &                      ECC, SEMI1, R1, G1
   11         FORMAT (/,' END QUINTUP2    T MT NM1 NMG NM3 E1 A1 R1 G1',
     &                                    F10.2,F6.2,3I6,F6.2,1P,3E10.2)
          END IF
      END IF
*
*       Include diagnostics for double triple ([[B,S],[B,S]]).
      IF (NAMEG(IMERGE).LT.0) THEN
          I1 = 2*IPAIR - 1
          J1 = 2*IPAIR
*       Obtain current merger index (note NAMEM(JM) = NAMEG(JG)).
          CALL FINDJ(J1,JJ,JG)
          IM = 1
          JM = 1
*       Determine original merger indices of the two inner binaries.
          DO 16 K = 1,NMERGE
              IF (NAMEM(K).EQ.NAME(I) + 2*NZERO) IM = K
              IF (NAMEM(K).EQ.NAMEG(JG)) JM = K
   16     CONTINUE
          AI = -0.5*(CM(1,IM) + CM(2,IM))/HM(IM)
          AJ = -0.5*(CM(1,JM) + CM(2,JM))/HM(JM)
          PMIN = SEMI1*(1.0 - ECC)
          ZM = BODY(I)*SMU
          if(rank.eq.0)
     &    WRITE (6,18)  TIME+TOFF, ZM, NAME(I1), NAMEG(IM), -NAMEM(JM),
     &                  NAMEG(JM), ECC, AI, AJ, R(IPAIR), SEMI1, PMIN,
     &                  G1
   18     FORMAT (/,' END HITRIP    T MT NM E1 AI AJ R1 A1 PM G1 ',
     &                              F9.2,F6.2,4I6,F6.2,1P,6E10.2)
      END IF
*
*       Check optional diagnostics for hierarchy.
      IF ((KZ(18).EQ.1.OR.KZ(18).EQ.3).AND.KSTAR(IMERGE).LE.20) THEN
          CALL HIARCH(IPAIR)
      END IF
*
*       Save neighbours for correction procedure.
      NNB = LIST(1,I) + 1
      DO 6 L = 2,NNB
          J = LIST(L,I)
          JPERT(L) = J
    6 CONTINUE
*
*       Ensure that c.m. coordinates are known to highest order.
      CALL XVPRED(I,0)
*
*       Predict neighbour coordinates & velocities (XDOT used by FPOLY1).
      DO 7 L = 2,NNB
          J = JPERT(L)
          CALL XVPRED(J,0)
    7 CONTINUE
*
*       Obtain current coordinates & velocities and specify KS components.
      CALL RESOLV(IPAIR,2)
      ICOMP = 2*IPAIR - 1
      JCOMP = ICOMP + 1
*
*       Initialize mass, coordinates and velocities for new c.m. body.
      BODY(I) = BODY(ICOMP)
      DO 8 K = 1,3
          X(K,I) = X(K,ICOMP)
          XDOT(K,I) = XDOT(K,ICOMP)
          X0DOT(K,I) = XDOT(K,ICOMP)
    8 CONTINUE
*
*       Add outer component to perturber list.
      JPERT(1) = JCOMP
*
*       Sum first part of potential energy correction due to tidal effect.
      JLIST(1) = ICOMP
      CALL NBPOT(1,NNB,POT1)
*
*       Find correct location of ghost particle using identification name.
      ICM = I
      JCOMP1 = JCOMP
      DO 10 I = 1,NTOT
          IF (BODY(I).EQ.0.0D0.AND.NAME(I).EQ.NAMEG(IMERGE)) JCOMP1 = I
   10 CONTINUE
*
*       Regularize two-body configuration if JCOMP1 cannot be identified.
      IF (JCOMP.EQ.JCOMP1) THEN
          if(rank.eq.0)
     &    WRITE (6,12)  IMERGE, NAMEG(IMERGE), JCOMP
   12     FORMAT (/,5X,'WARNING!    RESET2    JCOMP NOT IDENTIFIED ',
     &                       '   IM =',I3,'  NAMEG =',I6,'  JCOMP =',I6)
          GO TO 100
      END IF
*
*       Initialize basic variables for ghost and new c.m. (JCOMP -> JCOMP1).
      J1 = JCOMP1
      J = JCOMP
   13 T0(J1) = TIME
      BODY(J1) = BODY(J)
      DO 14 K = 1,3
          X(K,J1) = X(K,J)
          X0(K,J1) = X(K,J)
          XDOT(K,J1) = XDOT(K,J)
          X0DOT(K,J1) = XDOT(K,J)
   14 CONTINUE
      IF (J.EQ.JCOMP) THEN
          J1 = ICM
          J = ICOMP
          GO TO 13
      END IF
*
*       Restore masses, coordinates & velocities of inner binary.
      BODY(ICOMP) = CM(1,IMERGE)
      BODY(JCOMP) = CM(2,IMERGE)
      ZM = -BODY(ICOMP)/(BODY(ICOMP) + BODY(JCOMP))
*
*       Begin with second component since ICOMP holds new c.m. variables.
      I = JCOMP
      DO 20 KCOMP = 1,2
          DO 15 K = 1,3
              X(K,I) = X(K,ICOMP) + ZM*XREL(K,IMERGE)
              X0DOT(K,I) = X0DOT(K,ICOMP) + ZM*VREL(K,IMERGE)
              XDOT(K,I) = X0DOT(K,I)
*       Note that XDOT is only needed for improved polynomials of JCOMP.
   15     CONTINUE
          I = ICOMP
          ZM = BODY(JCOMP)/(BODY(ICOMP) + BODY(JCOMP))
   20 CONTINUE
*
*       Copy KS variables for inner binary (small TDOT2 near apo/peri).
      I1 = 2*IPAIR - 1
      T0(I1) = TIME
      LIST(1,I1) = 1
      H(IPAIR) = HM(IMERGE)
      R(IPAIR) = 0.0D0
      TDOT2(IPAIR) = 0.0D0
      DO 30 K = 1,4
          U(K,IPAIR) = UM(K,IMERGE)
          U0(K,IPAIR) = U(K,IPAIR)
          UDOT(K,IPAIR) = UMDOT(K,IMERGE)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
   30 CONTINUE
*
*       Add #JCOMP1 to all neighbour lists containing ICM.
      JLIST(1) = JCOMP1
      CALL NBREST(ICM,1,NNB)
*
*       Add #ICM to neighbour list of #JCOMP1.
      JLIST(1) = JCOMP1
      JPERT(1) = ICM
      CALL NBREST(ICM,1,1)
*
*       Form new neighbour list for #ICM (NB! NBREST(JCOMP1,1,1) skips).
      RS0 = RS(ICM)
      CALL NBLIST(ICM,RS0)
*
*       Get new list for #JCOMP1 (NNB = 1 from MERGE2/NBREM for ghost c.m.).
      RS0 = RS(JCOMP1)
      CALL NBLIST(JCOMP1,RS0)
*
*       Initialize force polynomial for outer component using resolved c.m.
      CALL FPOLY1(JCOMP1,JCOMP1,0)
      CALL FPOLY2(JCOMP1,JCOMP1,0)
*
*       Initialize c.m. polynomials and activate inner binary.
      CALL KSIN2(3)
*
*       Restore original name of inner hierarchy (c.m. NAME set in MERGE2).
      NAME(ICM) = NAME(ICM) + 2*NZERO
*
*       Locate position of inner binary in merger table (IM = 1 for safety).
      IM = 1
      DO 35 K = 1,NMERGE
          IF (NAMEM(K).EQ.NAME(ICM)) IM = K
   35 CONTINUE
*
*       Determine inclination between inner relative motion and outer orbit.
      RV = 0.0
      DO 40 K = 1,3
          XX(K,1) = XREL(K,IM)
          XX(K,2) = 0.0
          XX(K,3) = X(K,JCOMP)
          VV(K,1) = VREL(K,IM)
          VV(K,2) = 0.0
          VV(K,3) = XDOT(K,JCOMP)
          RV = RV + XREL(K,IM)*VREL(K,IM)
   40 CONTINUE
      CALL INCLIN(XX,VV,X(1,ICM),XDOT(1,ICM),ANGLE)
*
*       Evaluate stability parameter from current elements (minimum ECC).
      SEMI0 = -0.5*BODY(ICOMP)/HM(IM)
      SEMI = -0.5*BODY(ICM)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(ICM)*SEMI)
      ECC1 = SQRT(ECC2)
*       Form inner eccentricity from two-body elements.
      RIN = SQRT(XREL(1,IM)**2 + XREL(2,IM)**2 + XREL(3,IM)**2)
      ECC2 = (1.0 - RIN/SEMI0)**2 + RV**2/(BODY(ICOMP)*SEMI0)
      ECC = SQRT(ECC2)
*
*       Evaluate the general stability function.
      IF (ECC1.LT.1.0) THEN
          NST = NSTAB(SEMI0,SEMI,ECC,ECC1,ANGLE,CM(1,IMERGE),
     &                                    CM(2,IMERGE),BODY(JCOMP))
          IF (NST.EQ.0) THEN
              PCRIT = 0.98*SEMI*(1.0 - ECC1)
              PCR = stability(CM(1,IMERGE),CM(2,IMERGE),BODY(JCOMP),
     &                                          ECC,ECC1,ANGLE)*SEMI0
              ISTAB = ISTAB + 1
              IF (rank.eq.0.and.ISTAB.LT.50) THEN
                  WRITE (6,41)  ECC, ECC1, SEMI0, SEMI, PCRIT, PCR
   41             FORMAT (' STABLE2    E E1 A A1 PCR PC99 ',
     &                                 2F7.3,1P,4E10.2)
              END IF
          ELSE
              PCRIT = 1.02*SEMI*(1.0 - ECC1)
          END IF
      ELSE
          PCRIT = stability(CM(1,IMERGE),CM(2,IMERGE),BODY(JCOMP),
     &                                      ECC,ECC1,ANGLE)*SEMI0
      END IF
*
*       Set critical pericentre distance for stability check.
      R0(IPAIR) = PCRIT
*
*       Form new force polynomials for perturbers inside 10*R1 (skip J > N).
      NNB1 = LIST(1,ICM) + 1
      RS2 = (10.0*R1)**2
      DO 46 L = 2,NNB1
          J = LIST(L,ICM)
          RIJ2 = 0.0
          DO 42 K = 1,3
              RIJ2 = RIJ2 + (X(K,ICM) - X(K,J))**2
   42     CONTINUE
          IF (RIJ2.LT.RS2.AND.J.NE.JCOMP1.AND.J.LE.N) THEN
              T0(J) = TIME
              DO 45 K = 1,3
                  X0DOT(K,J) = XDOT(K,J)
   45         CONTINUE
              CALL FPOLY1(J,J,0)
              CALL FPOLY2(J,J,0)
          END IF
   46 CONTINUE
*
*       Rename perturber list for routine NBPOT.
      JPERT(1) = JCOMP1
*
*       Restore Roche stage indicator (any ghost c.m. is OK).
      KSTAR(ICM) = KSTARM(IMERGE)
*
*       See whether the outer component is a single or composite particle.
      POT3 = 0.0D0
      POT4 = 0.0D0
      IF (JCOMP1.LE.N) GO TO 50
*
*       Restore component masses for outer binary.
      JPAIR = JCOMP1 - N
      BODY(2*JPAIR-1) = CM(3,IMERGE)
      BODY(2*JPAIR) = CM(4,IMERGE)
*
*       Update look-up times & radii and check possible Roche condition.
      IF (KZ(19).GE.3) THEN
          IF (KSTAR(JCOMP1).GT.0.AND.KSTAR(JCOMP1).LE.10) THEN
              CALL TRFLOW(JPAIR,DTR)
              TEV(JCOMP1) = TIME + DTR
              TMDOT = MIN(TEV(JCOMP1),TMDOT)
              TMDOT = MIN(TEV(2*JPAIR),TMDOT)
          END IF
      END IF
*
*       Obtain coordinates & velocities of unperturbed binary components.
      CALL RESOLV(JPAIR,1)
*
*       Select new perturbers and initialize polynomials for KS motion.
      CALL KSLIST(JPAIR)
      CALL KSPOLY(JPAIR,1)
*
*       Apply tidal correction for outer binary perturbers.
      JLIST(1) = 2*JPAIR - 1
      JLIST(2) = 2*JPAIR
      CALL NBPOT(2,NNB,POT3)
      JLIST(1) = JCOMP1
      CALL NBPOT(1,NNB,POT4)
*
*       Update the merger energy.
      EB2 = BODY(2*JPAIR-1)*BODY(2*JPAIR)*H(JPAIR)/BODY(JCOMP1)
      EMERGE = EMERGE - EB2
*
      E2 = E1/EB2
      EB2 = EB2/BE(3)
      DP = POT4 - POT3
      IF (rank.eq.0.and.KZ(15).GT.1) THEN
          WRITE (6,48)  JPAIR, H(JPAIR), BODY(2*JPAIR-1),
     &                  BODY(2*JPAIR), E2, EB2, R(JPAIR), GAMMA(JPAIR),
     &                  DP
   48     FORMAT (' END OUTER MERGER',I5,'  H =',F7.1,'  M =',2F7.4,
     &            '  E1 =',F6.3,'  EB2 =',F6.3,'  RB2 =',1PE8.1,
     &            '  G2 =',E8.1,'  DP =',E8.1)
      END IF
*
*       Include interaction of body #ICOMP & JCOMP with perturbers.
   50 JLIST(1) = ICOMP
      JLIST(2) = JCOMP
      CALL NBPOT(2,NNB,POT2)
*
*       Form square of c.m. velocity correction due to tidal effects.
*     VI2 = X0DOT(1,ICM)**2 + X0DOT(2,ICM)**2 + X0DOT(3,ICM)**2
      DPHI = (POT2 - POT1) + (POT4 - POT3)
*     CORR = 1.0 + 2.0*DPHI/(BODY(ICM)*VI2)
*     IF (CORR.LE.0.0D0) CORR = 0.0
*
*       Adjust c.m. velocity by net tidal energy correction.
*     DO 60 K = 1,3
*         X0DOT(K,ICM) = SQRT(CORR)*X0DOT(K,ICM)
*  60 CONTINUE
*
*       Modify the merger energy to maintain conservation.
      EB = BODY(2*IPAIR-1)*BODY(2*IPAIR)*H(IPAIR)/BODY(ICM)
      EMERGE = EMERGE - EB + DPHI
*
      E1 = E1/EB
      EB = EB/BE(3)
      IF (rank.eq.0.and.KZ(15).GT.1) THEN
          WRITE (6,65)  IMERGE, TIME+TOFF, BODY(2*IPAIR-1),
     &                  BODY(2*IPAIR), R1, SEMI1, EB, E1,
     &                  GAMMA(IPAIR), G1, NNB-1
   65     FORMAT (' END MERGE2',I3,'  T =',F8.2,'  M =',2F7.4,
     &            '  R1 =',1PE8.1,'  A1 =',E8.1,'  EB =',0PF6.3,
     &            '  E1 =',F6.3,'  GB =',1PE8.1,'  G =',0PF6.3,
     &            '  NB =',I3)
      END IF
*
*       Check Roche look-up time.
      IF (KSTAR(ICM).GT.0.AND.KSTAR(ICM).LE.10) THEN
          K = ICM - N
          CALL TRFLOW(K,DTR)
          TEV(ICM) = MIN(TEV(ICM),TIME + DTR)
          TMDOT = MIN(TEV(ICM),TMDOT)
      END IF
*
*       Reduce merger counter and update tables (unless last or only pair).
   70 NMERGE = NMERGE - 1
      DO 80 L = IMERGE,NMERGE
          L1 = L + 1
          HM(L) = HM(L1)
          NAMEG(L) = NAMEG(L1)
          NAMEM(L) = NAMEM(L1)
          KSTARM(L) = KSTARM(L1)
          DO 74 K = 1,3
              XREL(K,L) = XREL(K,L1)
              VREL(K,L) = VREL(K,L1)
   74     CONTINUE
          DO 75 K = 1,4
              CM(K,L) = CM(K,L1)
              UM(K,L) = UM(K,L1)
              UMDOT(K,L) = UMDOT(K,L1)
   75     CONTINUE
   80 CONTINUE
*
*       Examine merger list for possible escapers (retain up to 3 levels).
      DO 90 L = 1,NMERGE
          DO 85 J = 1,NPAIRS
              IF (NAMEM(L).EQ.NAME(N+J).OR.
     &            NAMEM(L).EQ.NAME(N+J) + 2*NZERO.OR.
     &            NAMEM(L).EQ.NAME(N+J) + 4*NZERO) GO TO 90
   85     CONTINUE
*       Remove tables for any merger not identified.
          IMERGE = L
          GO TO 70
   90 CONTINUE
*
*       Set IPHASE < 0 for new time-step list in INTGRT.
      IPHASE = -1
*
  100 RETURN
*
      END
