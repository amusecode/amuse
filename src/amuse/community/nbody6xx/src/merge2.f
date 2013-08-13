      SUBROUTINE MERGE2
*
*
*       Merging of double hierarchy.
*       ----------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      CHARACTER*11  WHICH1
*
*
*       COMMON variables for merge procedure
*       ************************************
*
*       ------------------------------------------------------------------
*       CM      Component masses of merged binary (2nd binary in CM(3&4)).
*       HM      Binding energy per unit mass of merged binary.
*       KSTARM  Roche stage indicator (KSTAR for 1st c.m.; ghost is OK).
*       NAMEG   Name of the associated ghost component.
*       NAMEM   Name of new c.m. (< 0) for identification of merger index.
*       NMERG   Total number of mergers.
*       NMERGE  Current number of merged binaries (maximum is MMAX).
*       UM      Regularized coordinates of merged binary.
*       UMDOT   Regularized velocity of merged binary.
*       VREL    Relative velocity of merged binary components.
*       XREL    Relative coordinates of merged binary components.
*       ------------------------------------------------------------------
*
*
*       Ensure that the hierarchy is retained as primary KS pair.
      IF (NAME(N+KSPAIR).GT.0) THEN
          K = KSPAIR
          KSPAIR = JCOMP - N
          JCOMP = N + K
      END IF
*
*       Select deepest level as primary in case of two hierarchies.
      IF (NAME(N+KSPAIR).LT.0.AND.NAME(JCOMP).LT.0) THEN
*       Note possibility of a quartet (or even quintet) and triple.
          IF (NAME(JCOMP).LT.NAME(N+KSPAIR)) THEN
              K = KSPAIR
              KSPAIR = JCOMP - N
              JCOMP = N + K
          END IF
      END IF
*
*       Set pair index & c.m. of inner binary and save outer component.
      IPAIR = KSPAIR
      I = N + IPAIR
      JCOMP1 = JCOMP
      ICOMP1 = I
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      SEMI = -0.5*BODY(I)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
*
*       Produce diagnostics for [[B,S],B] quintuplet or sextuplet system.
      IF (NAME(JCOMP).GT.NZERO) THEN
          RI2 = 0.0
          DO 2 K = 1,3
              RI2 = RI2 + (X(K,I) - RDENS(K))**2
    2     CONTINUE
          RI = SQRT(RI2)
          EB = 0.5*BODY(I1)*BODY(I2)/SEMI
          EB = EB*FLOAT(N-NPAIRS)/ZKIN
          EB = MIN(EB,999.9D0)
          ZM = (BODY(I) + BODY(JCOMP))*SMU
          CALL FINDJ(I,JG,IM)
*       Check possible extension of next look-up time (skip GR case).
          IF (KZ(19).GE.3.AND.KZ(27).LT.3) THEN
              JLIST(1) = I1
              JLIST(2) = I2
              JLIST(3) = 2*(JCOMP - N) - 1
              JLIST(4) = JLIST(3) + 1
              JLIST(5) = JG
              CALL NEWTEV(5,IX)
          END IF
          WHICH1 = ' QUINTUPLET'
          IF (JG.GT.N) WHICH1 = ' SEXTUPLET '
          if(rank.eq.0)
     &    WRITE (6,5)  WHICH1, TIME+TOFF, ZM, NAME(I1), NAME(I2),
     &                 NAME(JCOMP), RI, ECC, EB, SEMI,PCRIT,GAMMA(IPAIR)
    5     FORMAT (/,' NEW',A11,'   T MT NM1 NM2 NM3 RI E0 EB0 A0 PC G0',
     &                             F10.2,F6.2,3I6,2F6.2,F6.1,1P,3E10.2)
      END IF
*
*       Check for [[B,S],S] quartet or [[B,B],S] quintuplet.
      IF (NAME(JCOMP).GT.0.AND.NAME(JCOMP).LE.NZERO) THEN
          CALL FINDJ(I,JG,IM)
          IF (KZ(19).GE.3.AND.KZ(27).LT.3) THEN
              JLIST(1) = I1
              JLIST(2) = I2
              JLIST(3) = JCOMP
              JLIST(4) = JG
              CALL NEWTEV(4,IX)
          END IF
          ZM = (BODY(I) + BODY(JCOMP))*SMU
          IF (rank.eq.0.and.JG.LE.N.AND.JCOMP.LE.N) THEN
              WRITE (6,7)  TIME+TOFF, ZM, NAME(I1), NAME(I2),
     &                     NAME(JCOMP), ECC, SEMI, PCRIT, GAMMA(IPAIR)
    7         FORMAT (/,' NEW QUARTET    T MT NM1 NMG NM3 E0 A0 PC G0 ',
     &                                   F9.2,F6.2,3I6,F6.2,1P,3E10.2)
          ELSE
              if(rank.eq.0)
     &        WRITE (6,8)  TIME+TOFF, ZM, NAME(I1), NAME(I2),
     &                     NAME(JCOMP), ECC, SEMI, PCRIT, GAMMA(IPAIR)
    8         FORMAT (/,' NEW QUINTUP2    T MT NM1 NMG NM3 E0 A0 PC G0',
     &                                    F10.2,F6.2,3I6,F6.2,1P,3E10.2)
          END IF
      END IF
*
*       Include diagnostics for double triple as [[B,S],[B,S]].
      IF (NAME(JCOMP).LT.0) THEN
          JPAIR = JCOMP - N
          J1 = 2*JPAIR - 1
          CALL FINDJ(I1,JI,IM)
          CALL FINDJ(J1,JJ,JM)
          IF (KZ(19).GE.3.AND.kZ(27).LT.3) THEN
              JLIST(1) = I1
              JLIST(2) = I2
              JLIST(3) = J1
              JLIST(4) = J1 + 1
              JLIST(5) = JI
              JLIST(6) = JJ
              CALL NEWTEV(6,IX)
          END IF
          AI = -0.5*(CM(1,IM) + CM(2,IM))/HM(IM)
          AJ = -0.5*(CM(1,JM) + CM(2,JM))/HM(JM)
          ZM = (BODY(I) + BODY(JCOMP))*SMU
          GX = MAX(GAMMA(IPAIR),GAMMA(JPAIR))
          if(rank.eq.0)
     &    WRITE (6,10)  TIME+TOFF, ZM, NAME(I1), NAME(I2), NAME(J1),
     &                  NAME(2*JPAIR), ECC, AI, AJ, R(IPAIR), R(JPAIR),
     &                  PCRIT, GX
   10     FORMAT (/,' NEW HITRIP    T MT NM E0 AI AJ RI RJ PC GX ',
     &                              F9.2,F6.2,4I6,F6.2,1P,6E10.2)
      END IF
*
*       Ensure correct coordinates & velocities.
      CALL RESOLV(IPAIR,3)
*
*       Check optional diagnostics for hierarchy.
      IF ((KZ(18).EQ.1.OR.KZ(18).EQ.3).AND.KSTAR(I).LE.20) THEN
          CALL HIARCH(IPAIR)
      END IF
*
*       Increase merger counter and set current merger index.
      NMERG = NMERG + 1
      NMERGE = NMERGE + 1
      IMERGE = NMERGE
*
*       Save component masses and evaluate reduced mass of inner binary.
      CM(1,IMERGE) = BODY(2*IPAIR-1)
      CM(2,IMERGE) = BODY(2*IPAIR)
      CM(3,IMERGE) = 0.0D0
      ZMU = BODY(2*IPAIR-1)*BODY(2*IPAIR)/BODY(N+IPAIR)
*
*       Set current energy of any second perturbed binary (from XVPRED).
      IF (JCOMP.GT.N) THEN 
          CALL XVPRED(JCOMP,-1)
          JPAIR = JCOMP - N
          IF (LIST(1,2*JPAIR-1).GT.0) H(JPAIR) = HT
*       Ensure that outer binary components are resolved.
          CALL KSRES(JPAIR,J1,J2,0.0D0)
*       Enforce non-zero membership for potential correction in NBPOT.
          LIST(1,2*JPAIR-1) = 1
      ELSE
          CALL XVPRED(JCOMP,-1)
      END IF
*
*       Save the neighbours for correction procedure.
      NNB = LIST(1,I)
      DO 20 L = 1,NNB
          J = LIST(L+1,I)
          JPERT(L) = J
   20 CONTINUE
*
*       Retain basic KS variables for explicit restart at merge termination.
      HM(IMERGE) = H(IPAIR)
      DO 25 K = 1,4
          UM(K,IMERGE) = U(K,IPAIR)
          UMDOT(K,IMERGE) = UDOT(K,IPAIR)
   25 CONTINUE
*
*       Save stellar type.
      KSTARM(IMERGE) = KSTAR(I)
*
*       Predict outer component to highest order if not on the block.
      IF (TIME - T0(JCOMP1).GT.0.0D0) THEN
          CALL XVPRED(JCOMP1,-1)
      END IF
*
*       Set temporary KS components for hierarchical binary.
      ICOMP = 2*IPAIR - 1
      JCOMP = ICOMP + 1
*
*       Obtain potential energy with respect to inner components.
      JLIST(1) = ICOMP
      JLIST(2) = JCOMP
      CALL NBPOT(2,NNB,POT1)
*
*       Save relative configuration.
      DO 30 K = 1,3
          XREL(K,IMERGE) = X(K,ICOMP) - X(K,JCOMP)
          VREL(K,IMERGE) = X0DOT(K,ICOMP) - X0DOT(K,JCOMP)
*       Initialize primary velocity of JCOMP1 (needed in KSIN2).
          X0DOT(K,JCOMP1) = XDOT(K,JCOMP1)
   30 CONTINUE
*
*       Include interaction of inner c.m. & neighbours before making new KS.
      ICOMP = ICOMP1
      JLIST(1) = ICOMP
      CALL NBPOT(1,NNB,POT2)
*
*       Copy mass of intruder to second KS component.
      BODY(JCOMP) = BODY(JCOMP1)
*
*       Replace name of #JCOMP with c.m. name (temporary for diagnostics).
      NAME2 = NAME(JCOMP)
      NAME(JCOMP) = NAME(JCOMP1)
*
*       Form new c.m. and initialize KS variables (JPERT safe first call!).
      JCOMP = JCOMP1
      CALL KSIN2(1)
*
*       See whether modifications due to second binary are needed.
      POT3 = 0.0D0
      POT4 = 0.0D0
      IF (JCOMP1.LE.N) GO TO 50
*
*       Initialize unperturbed ghost binary of outer component.
      T0(2*JPAIR-1) = 1.0E+06
      LIST(1,2*JPAIR-1) = 0
*
*       Apply tidal correction for outer binary perturbers.
      JLIST(1) = 2*JPAIR - 1
      JLIST(2) = 2*JPAIR
      CALL NBPOT(2,NNB,POT3)
      JLIST(1) = JCOMP1
      CALL NBPOT(1,NNB,POT4)
*
*       Update the merger energy to maintain conservation.
      EB1 = BODY(2*JPAIR-1)*BODY(2*JPAIR)*H(JPAIR)/BODY(JCOMP1)
      EMERGE = EMERGE + EB1
*
*       Save component masses and initialize ghost components.
      CM(3,IMERGE) = BODY(2*JPAIR-1)
      CM(4,IMERGE) = BODY(2*JPAIR)
      BODY(2*JPAIR-1) = 0.0D0
      BODY(2*JPAIR) = 0.0D0
      LIST(1,JCOMP) = 0
*
*       Remove ghost from all neighbour lists.
   50 JLIST(1) = JCOMP1
      ICM = N + KSPAIR
      CALL NBREM(ICM,1,NNB)
*
*       Specify JCOMP1 as ghost of zero mass.
      BODY(JCOMP1) = 0.0D0
*
*       Initialize integration variables to prevent spurious predictions.
      DO 60 K = 1,3
          X0DOT(K,JCOMP1) = 0.0D0
          XDOT(K,JCOMP1) = 0.0D0
          F(K,JCOMP1) = 0.0D0
          FDOT(K,JCOMP1) = 0.0D0
          D2(K,JCOMP1) = 0.0D0
          D3(K,JCOMP1) = 0.0D0
          D2R(K,JCOMP1) = 0.0D0
          D3R(K,JCOMP1) = 0.0D0
   60 CONTINUE
*
*       Set large value of T0 which avoids integration of ghost particle.
      T0(JCOMP1) = 1.0E+06
*       Set large X0 & X to avoid perturber selection (no escape removal).
      X0(1,JCOMP1) = 1.0E+06
      X(1,JCOMP1) = 1.0E+06
*
*       Initialize c.m. & KS polynomials.
      ICOMP = 2*IPAIR - 1
      JCOMP = ICOMP + 1
      CALL KSIN2(2)
*
*       Define large negative c.m. name for identification & termination.
      NAME(ICM) = NAME(ICM) - 3*NZERO
*
*       Set c.m. & ghost names for merger identification (include escape).
      NAMEM(IMERGE) = NAME(ICM)
      NAMEG(IMERGE) = NAME(JCOMP1)
      NAME(JCOMP) = NAME2
*
*       Copy stability limit for termination test A(1 - E) < R0 in KSINT.
      R0(IPAIR) = PCRIT 
*
*       Update merger energy to maintain conservation.
      DPHI = (POT2 - POT1) + (POT4 - POT3)
      EMERGE = EMERGE + ZMU*HM(IMERGE) + DPHI
*
*       Set IPHASE < 0 for new time-step list in routine INTGRT.
      IPHASE = -1
*
      RETURN
*
      END
