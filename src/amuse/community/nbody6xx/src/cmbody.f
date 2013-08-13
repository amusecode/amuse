      SUBROUTINE CMBODY(ENERGY,NSYS)
*
*
*       Formation of c.m. body by collision.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX4=4*NMX)
      COMMON/RCLOSE/  RIJ(4,4),RCOLL4,QPERI4,SIZE4(4),ECOLL4,IP(4)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/EBSAVE/  EBS
      REAL*8  CM(6),A0(3),A2(3)
      CHARACTER*8  WHICH1
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
*
*
*       Distinguish between chain and triple or quad case (ICH > 0 or = 0).
      IF (IPHASE.EQ.9) THEN
          ICH = 1
      ELSE
*       Activate collision indicator (otherwise done in CHTERM).
          ICH = 0
          IPHASE = 9
      END IF
*
*       Specify global indices of subsystem (membership: NSYS = 2 - 5).
      IF (NSYS.EQ.2) THEN
*
*       Define discrete time for prediction & new polynomials (T <= TBLOCK).
          I = N + KSPAIR
          DT = 0.1*STEP(I)
          IF (DT.GT.2.4E-11) THEN
              TIME2 = TIME - TPREV
              CALL STEPK(DT,DTN)
              TIME = TPREV + INT((TIME2 + DT)/DTN)*DTN
              TIME = MIN(TBLOCK,TIME)
          ELSE
              TIME = MIN(T0(I) + STEP(I),TBLOCK)
          END IF
          TIME0 = TIME
*
*       Check for hierarchical configuration.
          I1 = 2*KSPAIR - 1
          I2 = I1 + 1
          JCL = 0
          NP1 = LIST(1,I1) + 1
          SEMI = -0.5*BODY(I)/H(KSPAIR)
          RX = 100.0
          DO 5 L = 2,NP1
              J = LIST(L,I1)
              RIJ2 = 0.0
              VIJ2 = 0.0
              RDOT = 0.0
              A12 = 0.0
              A22 = 0.0
              A1A2 = 0.0
              DO 2 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
                  VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
                  RDOT = RDOT + (X(K,I) - X(K,J))*(XDOT(K,I) -XDOT(K,J))
                  K1 = K + 1
                  IF (K1.GT.3) K1 = 1
                  K2 = K1 + 1
                  IF (K2.GT.3) K2 = 1
                  A0(K) = (X(K1,I1)-X(K1,I2))*(XDOT(K2,I1)-XDOT(K2,I2))
     &                  - (X(K2,I1)-X(K2,I2))*(XDOT(K1,I1)-XDOT(K1,I2))
                  A2(K) = (X(K1,J) - X(K1,I))*(XDOT(K2,J) - XDOT(K2,I))
     &                  - (X(K2,J) - X(K2,I))*(XDOT(K1,J) - XDOT(K1,I))
                  A12 = A12 + A0(K)**2
                  A22 = A22 + A2(K)**2
                  A1A2 = A1A2 + A0(K)*A2(K)
    2         CONTINUE
              RIP = SQRT(RIJ2)
              A1 = 2.0/RIP - VIJ2/(BODY(I) + BODY(J))
              A1 = 1.0/A1
              ECC2 = (1.0 - RIP/A1)**2 +
     &                                  RDOT**2/(A1*(BODY(I) + BODY(J)))
              IF (1.0/A1.GT.0.5/RMIN) THEN
                  ECC1 = SQRT(ECC2)
                  RP1 = A1*(1.0 - ECC1)
                  ECC = 1.0 - R(KSPAIR)/SEMI
                  RA = SEMI*(1.0 + ECC)
                  SR = RP1/RA
                  GA = 2.0*BODY(J)*(RA/RP1)**3/BODY(I)
*       Determine inclination (8 bins of 22.5 degrees).
                  FAC = A1A2/SQRT(A12*A22)
                  ANGLE = 360.0*ACOS(FAC)/TWOPI
                  if(rank.eq.0)
     &            WRITE (6,4)  KSPAIR, NAME(J), H(KSPAIR), ECC, SEMI,
     &                         A1, RP1, GA, ECC1, SR, ANGLE
    4             FORMAT (' HIERARCHY:   KS NMJ H E A0 A1 RP GA E1 SR ',
     &                    'IN',2I6,F7.0,F9.5,1P,4E9.1,0P,F6.2,F6.1,F7.1)
              ELSE
                  IF (RIP.LT.RX) THEN
                      RX = RIP
                      JX = J
                      ECCX = SQRT(ECC2)
                      SEMIX = A1
                      RDX = RDOT/RX
                  END IF
              END IF
*
*       Select closest single body inside 0.5*RMIN as KS component.
              IF (RIP.LT.0.5*RMIN.AND.J.LE.N) THEN
                  IF (JCL.GT.0) THEN
                      IF (RIP.GT.RIP0) GO TO 5
                      JCL = J
                      RIP0 = RIP
                  ELSE
                      JCL = J
                      RIP0 = RIP
                  END IF
              END IF
    5     CONTINUE
*
          IF (rank.eq.0.and.RX.LT.20.0*RMIN) THEN
              WRITE (6,200)  NAME(JX), ECCX, RX, SEMIX*(1.0 - ECCX),
     &                       RDX, SEMI
  200         FORMAT (' PERTURBER:    NM E1 RX PM RD A0 ',
     &                                I6,F7.3,1P,4E9.1)
          END IF
*
*       Search for evidence of recent regularization.
          NAM1 = NAME(2*KSPAIR-1)
          NAM2 = NAME(2*KSPAIR)
          NNB = LISTD(1)
          DO 7 K = 2,NNB+1
              IF (LISTD(K).EQ.NAM1.OR.LISTD(K).EQ.NAM2) THEN
                  if(rank.eq.0)
     &            WRITE (6,6)  NAM1, NAM2, LISTD(K), K
    6             FORMAT (' KS REMNANT:    NAM LISTD K  ',3I6,I4)
              END IF
    7     CONTINUE
*
*       Predict body #JCL to current time in case of no collision.
          IF (JCL.GT.0) CALL XVPRED(JCL,-1)
*
*       Ensure orbit is at pericentre (perturbed hyperbolic case is OK).
          SEMI = -0.5*BODY(N+KSPAIR)/H(KSPAIR)
          IF (R(KSPAIR).GT.SEMI.AND.SEMI.GT.0.0) THEN
              CALL KSAPO(KSPAIR)
              CALL KSPERI(KSPAIR)
*       Restore quantized time to avoid small STEP (KSTERM needs T0 = TIME).
              TIME = TIME0
              T0(I1) = TIME
          END IF
*
*       Save collision distance and VINF before any common envelope stage.
          RCOLL = R(KSPAIR)
          VINF = 0.0
          IF (H(KSPAIR).GT.0.0) VINF = SQRT(2.0*H(KSPAIR))*VSTAR
          ECC = 1.0 - R(KSPAIR)/SEMI
*
*       Include special procedure for common envelope stage with mass loss.
          IF (KZ(19).GE.3) THEN
              K1 = KSTAR(I1)
              K2 = KSTAR(I2)
              ICASE = KTYPE(K1,K2)
              IF(ICASE.GT.100)THEN
                  IQCOLL = 4
                  CALL EXPEL(I1,I2,ICASE)
                  IF (ICASE.LT.0) GO TO 100
*       Treat collision as before in case of CE without coalescence.
                  ICOMP = I1
              END IF
          END IF
*
*       Update body #JCL to current time for new KS with combined c.m.
          IF (JCL.GT.0) THEN
*             CALL XVPRED(JCL,-1)
              T0(JCL) = TIME
              CALL DTCHCK(TIME,STEP(JCL),DTK(40))
              DO 8 K = 1,3
                  X0DOT(K,JCL) = XDOT(K,JCL)
                  X0(K,JCL) = X(K,JCL)
    8         CONTINUE
          END IF
*
*       Check diagnostics of degenerate binary (skip case of velocity kick).
          IF(KZ(8).GT.3.AND.MAX(KSTAR(I1),KSTAR(I2)).GE.10)THEN
              IF(KSTAR(I1).LE.12)THEN
                  CALL DEGEN(KSPAIR,KSPAIR,5)
              END IF
          END IF
*
*       Check optional diagnostics for final stage of binary evolution.
          IF (KZ(8).GT.3) THEN
              CALL BINEV(KSPAIR)
          END IF
*
*       Save binding energy (BODY(I2) = 0 is OK).
          EB = BODY(I1)*BODY(I2)*H(KSPAIR)/BODY(I)
          WHICH1 = ' BINARY '
          IF (H(KSPAIR).GT.0.0) THEN
              WHICH1 = ' HYPERB '
              NHYP = NHYP + 1
          END IF
*
*       Terminate KS pair and set relevant indices for collision treatment.
          KSTARI = KSTAR(I)
          T0(I1) = TIME
          SEMI = -0.5*BODY(I)/H(KSPAIR)
          TK = DAYS*SEMI*SQRT(SEMI/BODY(I))
          CALL DTCHCK(TIME,STEP(I1),DTK(40))
          CALL KSTERM
          I1 = 2*NPAIRS + 1
          I2 = I1 + 1
          I3 = 0
          ICOMP = I1
          DMIN2 = MIN(DMIN2,RCOLL)
      ELSE
*       Note JLIST(1->NCH) contains global indices (JLIST(4)=0 for NCH=3).
          I1 = JLIST(1)
          I2 = JLIST(2)
          I3 = JLIST(3)
          I4 = JLIST(4)
          IF (NSYS.GT.4) I5 = JLIST(5)
          IQCOLL = 5
          KSTARI = 0
          VINF = 0.0
          ECC = 1.0 + 2.0*EBS*DMINC/(BODY(I1)*BODY(I2))
          ECC = MAX(ECC,0.001D0)
          IF (EBS.GT.0) THEN
              HI = EBS*(BODY(I1) + BODY(I2))/(BODY(I1)*BODY(I2))
              VINF = SQRT(2.0*HI)*VSTAR
          END IF
*
*       Set new quantized time (note: suppress if problems in CHTERM).
          TIME = TBLOCK
*
*       Include special treatment for common envelope stage inside chain.
          IF (ICH.GT.0.AND.KZ(19).GE.3) THEN
              ZM1 = BODY(I1)*ZMBAR
              ZM2 = BODY(I2)*ZMBAR
              R1 = RADIUS(I1)*SU
              R2 = RADIUS(I2)*SU
              ECC = 1.0 + 2.0*EBS*DMINC/(BODY(I1)*BODY(I2))
              ECC = MAX(ECC,0.0D0)
              if(rank.eq.0)
     &        WRITE (86,9)  TPHYS, NAME(I1), NAME(I2), KSTAR(I1),
     &                      KSTAR(I2), ZM1, ZM2, R1, R2, DMINC*SU, ECC
    9         FORMAT (' COLL:    TPH NAM K* M R* QP E ',
     &                           F8.1,2I7,2I4,2F6.2,3F7.1,F9.5)
*
              K1 = KSTAR(I1)
              K2 = KSTAR(I2)
              ICASE = KTYPE(K1,K2)
              IF(ICASE.GT.100)THEN
                  IQCOLL = 6
                  CALL EXPEL2(I1,I2,ICASE)
*       Decide between chain restart, coalescence or collision.
                  IF (ICASE.GT.0) THEN
*       Adopt negative membership and reverse NSYS to denote chain restart.
                      NCH = -NSYS
                      NSYS = -NSYS
                      GO TO 100
                  END IF
*       Check for coalescence (one or two bodies of zero mass).
                  IF (BODY(I1).EQ.0.0D0.OR.BODY(I2).EQ.0.0D0) THEN
*       Include rare case of two massless stars (defined by KSTAR = 15).
                      IF (BODY(I1) + BODY(I2).EQ.0.0D0) THEN
                          RI = SQRT(X(1,I1)**2 + X(2,I1)**2 +X(3,I1)**2)
                          VI = SQRT(XDOT(1,I1)**2 + XDOT(2,I1)**2 +
     &                                              XDOT(3,I1)**2)
                          DTMAX = DTK(1)
                          CALL DTCHCK(TIME,DTMAX,DTK(40))
                          if(rank.eq.0)
     &                    WRITE (6,10)  KSTAR(I1), KSTAR(I2), DTMAX
   10                     FORMAT (' MASSLESS CM   K* DTX ',2I4,1P,E10.2)
                          GO TO 32
                      END IF
*       Reduce the membership (< 0 for SETSYS) and remove ghost from chain.
                      NCH = -(NSYS - 1)
                      JLIST(1) = I1
                      IF (BODY(I1).EQ.0.0D0) JLIST(1) = I2
                      DO 11 L = 2,NSYS-1
                          JLIST(L) = JLIST(L+1)
   11                 CONTINUE
                      GO TO 100
                  END IF
*       Treat collision as before in case of CE without coalescence.
                  ICOMP = I1
                  JLIST(1) = I1
                  JLIST(2) = I2
              END IF
          END IF
      END IF
*
*       Obtain mass loss and evolution epoch of composite star.
      DM = 0.0D0
      IF (KZ(19).GE.3) THEN
          CALL MIX(I1,I2,DM)
          ICOMP = I1
*       Note possible switching of I1 and I2 (cf. JLIST).
      END IF
*
*       Define global c.m. coordinates & velocities from body #I1 & I2.
      ZM = BODY(I1) + BODY(I2)
      DO 12 K = 1,3
          CM(K) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/ZM
          CM(K+3) = (BODY(I1)*XDOT(K,I1) + BODY(I2)*XDOT(K,I2))/ZM
   12 CONTINUE
*
*	Set T0 = TIME for correct potential energy correction in FCORR.
      IF (ICH.GT.0) THEN
          DO 14 L = 1,NCH
              J = JLIST(L)
              T0(J) = TIME
              CALL DTCHCK(TIME,STEP(J),DTK(40))
   14     CONTINUE
      END IF
*
*       Ensure the heaviest body is new progenitor (only ICOMP is needed).
      IF (BODY(I2).GT.BODY(I1)) THEN
          I1S = I1
          I1 = I2
          I2 = I1S
          ICOMP = I1
          JCOMP = I2
      END IF
*
*       Evaluate potential energy with respect to colliding bodies.
      IF (NSYS.EQ.2) THEN
*       Copy perturber list to JPERT.
          NNB = LIST(1,I1)
          DO 15 L = 1,NNB
              JPERT(L) = LIST(L+1,I1)
   15     CONTINUE
          JLIST(1) = I1
          JLIST(2) = I2
*       Replace second old KS component temporarily by arbitrary body.
          JPERT(1) = N
          IF (I2.EQ.N) JPERT(1) = N + 1
          CALL NBPOT(2,NNB,POT1)
      ELSE
*       Obtain differential effect on #I1 & #I2 due to other members.
          DO 16 L = 3,NSYS
              JPERT(L-2) = JLIST(L)
   16     CONTINUE
          NP = NSYS - 2
          CALL NBPOT(2,NP,POT1)
*       Compensate EGRAV for any chain mass loss.
          IF (DM.GT.0.0D0) THEN
              EG0 = EGRAV
              DO 18 L = 1,NSYS
                  J = JLIST(L)
                  IF (J.EQ.I1.OR.J.EQ.I2) GO TO 18
                  RIJ2 = 0.0
                  DO 17 K = 1,3
                      RIJ2 = RIJ2 + (CM(K) - X(K,J))**2
   17             CONTINUE
                  EGRAV = EGRAV - DM*BODY(J)/SQRT(RIJ2)
   18         CONTINUE
              if(rank.eq.0)
     &        WRITE (6,19)  NSYS, EGRAV, EGRAV - EG0
   19         FORMAT (' CMAIN MASS LOSS    NSYS EGRAV DEGR ',
     &                                     I4,1P,2E10.2)
          END IF
      END IF
*
*       Create new body from c.m. and initialize zero mass ghost in #I2.
      ZM1 = BODY(I1)*ZMBAR
      ZM2 = BODY(I2)*ZMBAR
      BODY(I1) = ZM
      BODY(I2) = 0.0D0
      SPIN(I1) = (SPIN(I1) + SPIN(I2))*(1.0 - DM/ZM)
      T0(I2) = TADJ + DTADJ 
      IF (KZ(23).EQ.0.OR.RTIDE.GT.1000.0*RSCALE) T0(I2) = 1.0D+10
      DTMAX = DTK(1)
      CALL DTCHCK(TIME,DTMAX,DTK(40))
      STEP(I2) = DTMAX
      RI = SQRT((X(1,I2) - RDENS(1))**2 + (X(2,I2) - RDENS(2))**2
     &                                  + (X(3,I2) - RDENS(3))**2)
      VI = SQRT(XDOT(1,I2)**2 + XDOT(2,I2)**2 + XDOT(3,I2)**2)
      NAME1 = NAME(I1)
      NAME2 = NAME(I2)
*
      DO 20 K = 1,3
          X(K,I1) = CM(K)
          X0(K,I1) = CM(K)
          XDOT(K,I1) = CM(K+3)
          X0DOT(K,I1) = CM(K+3)
*       Ensure that ghost will escape next output (far from fast escapers).
          X0(K,I2) = 1000.0*RSCALE*(X(K,I2) - RDENS(K))/RI
          X(K,I2) = X0(K,I2)
          X0DOT(K,I2) = SQRT(0.004*ZMASS/RSCALE)*XDOT(K,I2)/VI
          XDOT(K,I2) = X0DOT(K,I2)
          F(K,I2) = 0.0D0
          FDOT(K,I2) = 0.0D0
          D2(K,I2) = 0.0D0
          D3(K,I2) = 0.0D0
   20 CONTINUE
*
*       Set appropriate parameters for coalescing GR binaries.
      IF (KSTAR(I1).GE.13.AND.KZ(28).GT.0) THEN
          TEV(I1) = 1.0D+10
      END IF
      IF (KSTAR(I2).GE.13.AND.KZ(28).GT.0) THEN
          TEV(I2) = 1.0D+10
          KSTAR(I2) = 0
      END IF
*
*       Refresh index of new dominant body in case of switch in routine MIX.
      JLIST(1) = I1
*       Obtain potential energy w.r.t. new c.m. and apply tidal correction.
      IF (NSYS.EQ.2) THEN
          CALL NBPOT(1,NNB,POT2)
      ELSE
          CALL NBPOT(1,NP,POT2)
      END IF
      DP = POT2 - POT1
      ECOLL = ECOLL + DP
*
*       Remove the ghost particle from perturber lists containing #I1.
      IF (NSYS.EQ.2) THEN
          JPERT(1) = I2
          JLIST(1) = I2
          CALL NBREM(I1,1,NNB)
*       Remove ghost from list of I1 (use NTOT as dummy here).
          JPERT(1) = I1
          CALL NBREM(NTOT,1,1)
      ELSE
*       Determine index and set neighbour membership of original chain c.m.
          JLIST(2) = I2
          NNB = 0
          DO 21 L = 1,NSYS
              J = JLIST(L)
              IF (LIST(1,J).GT.NNB) THEN
                  NNB = LIST(1,J)
                  ICM = J
              END IF
   21     CONTINUE
*
*       Update neighbour lists of current chain c.m. and remove ghost I2.
          DO 22 L = 1,NNB
              JPERT(L) = LIST(L+1,ICM)
   22     CONTINUE
          CALL NBREST(ICM,NSYS,NNB)
          JLIST(1) = I2
          JPERT(1) = ICM
          CALL NBREM(NTOT,1,1)
      END IF
*
*       Include correction procedure in case of mass loss (cf routine MIX).
      IF (KZ(19).GE.3.AND.DM.GT.0.0D0) THEN
*
*       Determine neighbour list for chain member.
          IF (NSYS.GT.2) THEN
              RSI = RSCALE*(10.0/FLOAT(N - NPAIRS))**0.3333
              CALL NBLIST(I1,RSI)
          END IF
*
          NNB = LIST(1,I1)
          ILIST(1) = NNB
          DO 24 L = 2,NNB+1
              ILIST(L) = LIST(L,I1)
   24     CONTINUE
*
*       Reduce mass of composite body and update total mass (check SN mass).
          BODY(I1) = ZM - DM
          BODY(I1) = MAX(BODY(I1),0.0D0)
          ZMASS = ZMASS - DM
          ZM = ZM - DM
*
*       Perform total force & energy corrections (new polynomial set later).
          KW = KSTAR(I1)
          CALL FCORR(I1,DM,KW)
*
*       Initialize new polynomials of neighbours & #I for DM > 0.1 DMSUN.
          IF (DM*ZMBAR.GT.0.1) THEN
*
*       Include body #I at the end (counting from location #2; not KS case).
              NNB2 = NNB + 2
              ILIST(NNB2) = I1
              IF (NSYS.EQ.2) NNB2 = NNB2 - 1
*
*       Obtain new F & FDOT and time-steps.
              DO 30 L = 2,NNB2
                  J = ILIST(L)
                  DO 25 K = 1,3
                      X0DOT(K,J) = XDOT(K,J)
   25             CONTINUE
                  CALL FPOLY1(J,J,0)
                  CALL FPOLY2(J,J,0)
   30         CONTINUE
          END IF
          TPREV = TIME - STEPX
      END IF
*
*       Check creation of ghost(s) after collision of two white dwarfs.
   32 IF (KSTAR(I1).EQ.15) THEN
          T0(I1) = TADJ + DTADJ 
          STEP(I1) = DTMAX
          DO 35 K = 1,3
              X0(K,I1) = 1000.0*RSCALE*X(K,I1)/RI
              X(K,I1) = X0(K,I1)
              X0DOT(K,I1) = SQRT(0.004*ZMASS/RSCALE)*XDOT(K,I2)/VI
              XDOT(K,I1) = X0DOT(K,I1)
              F(K,I1) = 0.0D0
              FDOT(K,I1) = 0.0D0
              D2(K,I1) = 0.0D0
              D3(K,I1) = 0.0D0
   35     CONTINUE
*       Include case of two ghost stars.
          IF (BODY(I2).EQ.0.0D0.AND.I1.NE.I2) THEN
              I1 = I2
              GO TO 32
          END IF
*       Perform initialization of single particle in three-body chain.
          IF (NSYS.EQ.3) THEN
              ICOMP = I3
*       Make quick exit from routine CHAIN on zero membership.
              NSYS = 0
              NCH = 0
              ECH = 0.0
              GO TO 40
          END IF
*       Include treatment for larger chain here when needed (rare case).
          IF (NSYS.GT.3) THEN
              if(rank.eq.0)
     &        WRITE (6,38)  NAME(I1)
   38         FORMAT (' DANGER!    ZERO MASS IN CHAIN    NAME',I7)
              STOP
          END IF
      END IF
*
*       Decide appropriate path for each case.
      IF (NSYS.EQ.2) GO TO 40
      IF (NSYS.EQ.3) GO TO 45
*
*       Switch KS components if body #I3 & I4 is closer than #I1 & I3.
      IF (JLIST(7).LT.0) THEN
          I4 = I1
          I1S = I1
          I1 = JLIST(4)
          JLIST(4) = I1S
      END IF
      ICOMP = I4
*
*       Check KS case for new regularization with close hierarchical body.
   40 IF (NSYS.EQ.2) THEN
          IF (JCL.GT.0) THEN
              ICOMP = I1
              JCOMP = JCL
              CALL KSREG
              GO TO 80
          END IF
      END IF
*
*       Initialize force polynomial for new single, third or fourth body.
      CALL FPOLY1(ICOMP,ICOMP,0)
      CALL FPOLY2(ICOMP,ICOMP,0)
      IF (NSYS.EQ.5) THEN
          CALL FPOLY1(I5,I5,0)
          CALL FPOLY2(I5,I5,0)
          if(rank.eq.0)
     &    WRITE (6,42)  NAME(I5), STEP(I5)
   42     FORMAT (' 5-CHAIN!    NM DT ',I7,1P,E10.2)
      END IF
      IF (NSYS.EQ.0) GO TO 95
      IF (NSYS.LE.2) GO TO 80
*
*       Add kinetic energy from last body and check DMIN in TRIPLE or QUAD.
   45 IF (ICH.GT.0) THEN
          WHICH1 = '  CHAIN '
*       Specify new membership (< 0 for SETSYS) and remove ghost from chain.
          NCH = -(NSYS - 1)
          JLIST(1) = I1
          DO 50 L = 2,NSYS-1
              JLIST(L) = JLIST(L+1)
   50     CONTINUE
*       Copy well defined binding energy and skip explicit evaluation.
          EB = EBS
          CHCOLL = CHCOLL + EB
          GO TO 80
      ELSE IF (NSYS.EQ.3) THEN
          I = I3
          WHICH1 = ' TRIPLE '
          DMIN3 = MIN(DMIN3,RCOLL4)
          RCOLL = RCOLL4
          EB = EBS
      ELSE
          I = I4
          WHICH1 = '   QUAD '
          DMIN4 = MIN(DMIN4,RCOLL4)
          RCOLL = RCOLL4
          EB = EBS
      END IF
*
*       Set global components for new KS regularization (ICOMP < JCOMP).
      ICOMP = MIN(I1,I3)
      JCOMP = MAX(I1,I3)
*
*       Initialize new KS pair.
      CALL KSREG
*
*       Update energy loss & collision counters.
   80 ECOLL = ECOLL + EB
      E(10) = E(10) + EB + DP
      EGRAV = EGRAV + EB
*
*       Open the second coalescence unit #26 first time.
      IF (FIRST.AND.(IQCOLL.EQ.3.OR.KSTARI.GE.10)) THEN
          OPEN (UNIT=26,STATUS='UNKNOWN',FORM='FORMATTED',FILE='COAL2')
          FIRST = .FALSE.
*
*       Print cluster scaling parameters at start of the run.
          if(rank.eq.0)then
          WRITE (26,82)  RBAR, BODYM*ZMBAR, BODY1*ZMBAR, TSCALE,
     &                   NBIN0, NZERO
   82     FORMAT (/,4X,'MODEL:    RBAR =',F5.1,'  <M> =',F6.2,
     &                 '  M1 =',F6.1,'  TSCALE =',F6.2,
     &                 '  NB =',I4,'  N0 =',I6,//)
*
          WRITE (26,84)
   84     FORMAT ('    TIME  NAME  NAME  K1  K2  IQ  M1   M2',
     &            '   DM    R1     R2    r/Rc   R     ECC      P',/)
          end if
      END IF
*
*       Distinguish case of contact binary (i.e. coalescence).
      IF (IQCOLL.EQ.3.OR.KSTARI.GE.10) THEN
          NPOP(8) = NPOP(8) + 1
          NCOAL = NCOAL + 1
          if(rank.eq.0)then
          WRITE (6,85)  IQCOLL, NAME1, NAME2, ZM*SMU, RCOLL, EB, DP, ECC
   85     FORMAT (/,' BINARY COAL    IQCOLL =',I3,'  NAME =',2I6,
     &             '  M =',F6.2,'  RCOLL =',1P,E8.1,' EB =',E9.1,
     &             '  DP =',E9.1,'  E =',0P,F8.4)
*
          WRITE (26,86)  TTOT, NAME1, NAME2, KSTAR(I1), KSTAR(I2),
     &                   IQCOLL, ZM1, ZM2, DM*ZMBAR, R1, R2, RI/RC,
     &                   RCOLL*SU, ECC, TK
   86     FORMAT (1X,F7.1,2I6,3I4,3F5.1,2F7.2,F6.2,F7.2,F9.5,1P,E9.1)
          CALL FLUSH(26)
          end if
          GO TO 95
      END IF
*
      NPOP(8) = NPOP(8) + 1
      NCOLL = NCOLL + 1
*
      if(rank.eq.0)
     &WRITE (6,90)  WHICH1, NSYS, NAME1, NAME2, ZM*SMU, RCOLL, EB,
     &              VINF, ECC, DP
   90 FORMAT (/,A8,'COLLISION    NSYS =',I3,'  NAME =',2I6,
     &             '  M =',F6.2,'  RCOLL =',1P,E8.1,'  EB =',E9.1,
     &             '  VINF =',0P,F5.1,'  ECC =',F9.5,'  DP =',1P,E9.1)
*
*       Specify IPHASE < 0 for new sorting.
   95 IPHASE = -1
*
*       Reduce NSUB for chain (temporary increase by CHINIT before CHTERM).
  100 IF (ICH.GT.0) THEN
          NSUB = NSUB - 1
      END IF
*       Skip NSUB reduction for continuation of CHAIN (bug fix 26/8/06).
      TTOT = TIME + TOFF
*
      RETURN
*
      END
