      SUBROUTINE KSTIDE(IPAIR,QPERI)
*
*
*       Tidal or GR interaction of KS pair.
*       -----------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  DE(2)
      CHARACTER*8  WHICH1
      INTEGER  IS(2)
      SAVE  TTIDE,IONE
      DATA  ECCM,ECCM2,TTIDE,IONE  /0.002,0.00000399,0.0D0,0/
      SAVE SUM
      DATA SUM /0.0D0/
*
*
*       Define indices of KS components.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
*
*       Set c.m. index & reduced mass and increase event counter.
      I = N + IPAIR
      ZMU = BODY(I1)*BODY(I2)/BODY(I)
      NDISS = NDISS + 1
      RKS = R(IPAIR)
*
*       Form current semi-major axis & eccentricity (not at peri).
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(SEMI*BODY(I))
      ECC = SQRT(ECC2)
      AM0 = SEMI*(1.0D0 - ECC**2)
      PERI = SEMI*(1.0D0 - ECC)
      HI = H(IPAIR)
*
*       Distinguish between sequential circularization, standard chaos & GR.
      R1 = MAX(RADIUS(I1),RADIUS(I2))
      IF (KZ(27).EQ.1) THEN
          ZF = 4.0
          IF (ECC.GT.0.95) ZF = 50.0
*       Skip if penetration is not significant (prevents frequent calls).
          IF (ABS(QPERI - ZF*R1).LT.0.01*QPERI) GO TO 50
      ELSE IF (KZ(27).EQ.2) THEN
          IS(1) = KSTAR(I1)
          IS(2) = KSTAR(I2)
*       Obtain kinetic energy loss due to tidal interaction (DE > 0).
*        CALL TIDES(QPERI,BODY(I1),BODY(I2),RADIUS(I1),RADIUS(I2),IS,DE)
*
*       Evaluate chaos parameters and check for disruption.
          CALL CHAOS(IPAIR,I1,I2,QPERI,ECC,IS,ZMU,RKS,SEMI1,ECC1,IDIS)
          IF (IDIS.EQ.1) GO TO 45
*
*       Update KS variables for CHAOS or TERMINATED CHAOS (H > 0).
          IF (KSTAR(I).EQ.-1.OR.H(IPAIR).GT.0.0) THEN
*       Ensure rectification to determine true pericentre (bug fix 07/08).
              CALL KSRECT(IPAIR)
              CALL KSPERI(IPAIR)
              GO TO 10
          END IF
*       Skip on WIDE CHAOS or NEW SPIRAL.
          GO TO 45
      ELSE IF (KZ(27).EQ.3) THEN
          CALL TIDES3(QPERI,BODY(I1),BODY(I2),VSTAR,H(IPAIR),ECC,DE)
      END IF
*
*       Restore circularization index if needed (exit from CHAIN).
      IF (ECC.LE.ECCM.AND.KSTAR(I).LT.10) THEN
          KSTAR(I) = 10
          GO TO 50
      END IF
*
*       Consider sequential circularization or GR evolution.
      IF (KZ(27).EQ.1) THEN
*       Suppress the old PT procedure (DH => ECC from AM0 = const).
*         ECC2 = ECC**2 + 2.0D0*AM0*DH/BODY(I)
*         ECC2 = MAX(ECC2,ECCM2)
*       Accept circularized orbit if ACIRC < 4*R1 (use maximum radius).
          AM0 = SEMI*(1.0 - ECC**2)
          ECC1 = SQRT(ECCM2)
          ACIRC = AM0/(1.0 - ECCM2)
          IF (ACIRC.LT.ZF*R1) THEN
              SEMI1 = ACIRC
          ELSE
*       Obtain E1 from A1*(1 - E1**2) = AM0 using A1*(1 - E1) = 4*R1.
              ECC1 = AM0/(ZF*R1) - 1.0
              ECC1 = MAX(ECC1,ECCM)
*       Set new semi-major axis from angular momentum conservation.
              SEMI1 = AM0/(1.0 - ECC1**2)
          END IF
*       Form the corresponding energy change.
          DH = 0.5*BODY(I)*(1.0/SEMI - 1.0/SEMI1)
          DE(1) = -ZMU*DH
          DE(2) = 0.0
      ELSE
*       Include safety check on energy loss to prevent new SEMI < R.
          DH = -(DE(1) + DE(2))/ZMU
          IF (H(IPAIR) + DH.LT.-0.5*BODY(I)/R(IPAIR)) THEN
              DH = -0.5*BODY(I)/R(IPAIR) - H(IPAIR)
              DE(1) = -ZMU*DH
              DE(2) = 0.0
          END IF
          SEMI1 = -0.5*BODY(I)/(H(IPAIR) + DH)
          ECC1 = 1.0 - PERI/SEMI1
          ECC1 = MAX(ECC1,ECCM)
      END IF
*
*       Skip possible hyperbolic case.
      IF (H(IPAIR) + DH.GT.0.0) GO TO 50
*
*       Update total energy loss (and H after obtaining peri).
      ECOLL = ECOLL + (DE(1) + DE(2))
      E(10) = E(10) + (DE(1) + DE(2))
      EGRAV = EGRAV + (DE(1) + DE(2))
*       Determine pericentre variables U & UDOT by backwards integration.
      CALL KSPERI(IPAIR)
      H(IPAIR) = H(IPAIR) + DH
*
*       Print first and last energy change and update indicator.
      IF (KZ(27).EQ.1.AND.(KSTAR(I).EQ.0.OR.KSTAR(I).EQ.9)) THEN
          P = DAYS*SEMI1*SQRT(SEMI1/BODY(I))
          IF (KSTAR(I).EQ.0.AND.ECC1.GT.ECCM) THEN
              WHICH1 = 'NEW CIRC'
              KSTAR(I) = 9
              if(rank.eq.0)
     &        WRITE (6,8)  WHICH1, NAME(I1), NAME(I2), KSTAR(I1),
     &                     KSTAR(I2), TTOT, ECC, ECC1, P, SEMI1, R1
          ELSE IF (ECC1.LE.ECCM) THEN
              WHICH1 = 'END CIRC'
              KSTAR(I) = 10
              NCIRC = NCIRC + 1
              if(rank.eq.0)
     &        WRITE (6,8)  WHICH1, NAME(I1), NAME(I2), KSTAR(I1),
     &                     KSTAR(I2), TTOT, ECC, ECC1, P, SEMI1, R1
          END IF
    8     FORMAT (' ',A8,'    NAM K* T E0 EF P AF R* ',
     &                        2I6,2I4,F9.2,2F8.3,F7.1,1P,2E10.2)
      END IF
*
*       Set new pericentre.
   10 PERI1 = SEMI1*(1.0D0 - ECC1)
*
*       Form KS coordinate scaling factor from pericentre ratio.
      C1 = SQRT(PERI1/PERI)
*
*       Specify KS velocity scaling (conserved J, chaos or GR treatment).
      IF (KZ(27).EQ.1) THEN
          C2 = 1.0/C1
      ELSE
*       Note that PERI should not change in GR case (hence same C2).
          C2 = SQRT((BODY(I) + H(IPAIR)*PERI1)/(BODY(I) + HI*PERI))
      END IF
*
*       See whether circular orbit condition applies.
*     IF (ECC1.LE.ECCM) THEN
*         AM = SEMI1*(1.0D0 - ECC1**2)
*         C2 = SQRT(AM/AM0)/C1
*     END IF
*
*       Transform KS variables to yield the prescribed elements.
      R(IPAIR) = 0.0D0
      DO 15 K = 1,4
          U(K,IPAIR) = C1*U(K,IPAIR)
          UDOT(K,IPAIR) = C2*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
   15 CONTINUE
*
*       Form new perturber list after significant energy loss.
      NP0 = LIST(1,I1)
      IF (ABS(SEMI1/SEMI).LT.0.5) THEN
          CALL KSLIST(IPAIR)
*       Ensure that single perturber differs from #I itself.
          IF (NP0.EQ.0.AND.LIST(1,I1).GT.0) THEN
              IF (LIST(2,I1).EQ.I) THEN
                  LIST(2,I1) = IFIRST
              END IF
          END IF
      END IF
*
*       Ensure chaos treatment at every pericentre by perturbed motion.
      IF (KZ(27).EQ.2.AND.LIST(1,I1).EQ.0) THEN
          LIST(1,I1) = 1
          LIST(2,I1) = N
          NP0 = 1
      END IF
*
*       Copy perturber list to JPERT and predict perturbers.
      NNB = LIST(1,I1)
      DO 16 L = 1,NNB
          JPERT(L) = LIST(L+1,I1)
          J = LIST(L+1,I1)
          CALL XVPRED(J,0)
   16 CONTINUE
*
*       Evaluate potential energy of perturbers.
      JLIST(1) = I1
      JLIST(2) = I2
      CALL NBPOT(2,NNB,POT1)
*
*       Re-initialize KS polynomials at pericentre for perturbed case.
      T0(I1) = TIME
      IF (NP0.GT.0.OR.ECC.LE.ECCM) THEN
          CALL RESOLV(IPAIR,1)
          IMOD = KSLOW(IPAIR)
          CALL KSPOLY(IPAIR,IMOD)
      END IF
*
*       Obtain potential energy wrt new binary and apply tidal correction.
      CALL NBPOT(2,NNB,POT2)
      DP = POT2 - POT1
      ECOLL = ECOLL + DP
*
*       Produce diagnostic output for interesting new case (first time).
      IF (ECC.GT.0.99.AND.ABS(ECC - ECC1).GT.0.01.AND.IONE.EQ.0) THEN
          if(rank.eq.0)
     &    WRITE (6,20)  NAME(I1), NAME(I2), SEMI1, ECC, ECC1, HI,
     &                  QPERI, DH, DP
   20     FORMAT (' NEW KSTIDE    NAM AF E0 EF HI QP DH DP ',
     &                            2I5,1P,E10.2,0P,2F8.3,F9.1,1P,3E10.2)
          TTIDE = TIME + TOFF
          IONE = IONE + 1
      END IF
      IF (TIME + TOFF.GT.TTIDE + DTADJ) IONE = 0
*
*       Check for hierarchical configuration with eccentric inner binary.
      IF (KZ(27).EQ.2.AND.ECC.GT.0.95.AND.HI.LT.0.0) THEN
          NP1 = LIST(1,I1) + 1
          DO 30 L = 2,NP1
              J = LIST(L,I1)
              RIJ2 = 0.0
              VIJ2 = 0.0
              RDOT = 0.0
              DO 25 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
                  VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
                  RDOT = RDOT + (X(K,I) - X(K,J))*(XDOT(K,I) -XDOT(K,J))
   25         CONTINUE
              RIP = SQRT(RIJ2)
              A1 = 2.0/RIP - VIJ2/(BODY(I) + BODY(J))
              A1 = 1.0/A1
              IF (1.0/A1.GT.0.1/SEMI) THEN
                  ECC2 = (1.0 - RIP/A1)**2 +
     &                                  RDOT**2/(A1*(BODY(I) + BODY(J)))
                  RP = A1*(1.0 - SQRT(ECC2))
                  RA = SEMI*(1.0 + ECC)
                  SR = RP/RA
                  ICIRC = -1
                  JCOMP = J
                  CALL INDUCE(IPAIR,EMAX,EMIN,ICIRC,TC,ANGLE,TG,EDAV)
                  if(rank.eq.0)
     &            WRITE (6,28)  NAME(J), H(IPAIR), SEMI, A1, RP, EDAV,
     &                          SQRT(ECC2), EMAX, SR
   28             FORMAT (' HIERARCHY    NMJ H A0 A1 RP EDAV E1 EX SR ',
     &                                   I7,F7.0,1P,4E9.1,0P,2F8.4,F6.1)
              END IF
   30     CONTINUE
      END IF
*
*       Set one unperturbed period for small apocentre perturbation (#27=1).
      GA = GAMMA(IPAIR)*(SEMI1*(1.0 + ECC1)/R(IPAIR))**3
      IF (GA.LT.GMIN.AND.KZ(27).EQ.1) THEN
          STEP(I1) = TWOPI*SEMI1*SQRT(SEMI1/BODY(I))
          LIST(1,I1) = 0
      END IF
*
*       Ensure T'' = 0 for pericentre test in KSINT & dissipation in UNPERT.
      IF (TDOT2(IPAIR).LT.0.0D0) THEN
          TDOT2(IPAIR) = 0.0D0
      END IF
*
*       Count any hyperbolic captures.
      IF (SEMI.LT.0.0.AND.SEMI1.GT.0.0) THEN
          NTIDE = NTIDE + 1
          QPS = QPERI/R1
          if(rank.eq.0)
     &    WRITE (6,35)  TIME+TOFF, NAME(I1), NAME(I2), ECC, ECC1, QPS,
     &                  SEMI1
   35     FORMAT (' NEW CAPTURE    T NM E EF QP/R* A1 ',
     &                             F9.2,2I6,2F9.4,1P,2E10.2)
      END IF
*
*       Record diagnostics for new synchronous orbit and activate indicator.
      IF (ECC.GT.ECCM.AND.ECC1.LT.ECCM.AND.KZ(27).LE.2) THEN
          NSYNC = NSYNC + 1
          ESYNC = ESYNC + ZMU*H(IPAIR)
          KSTAR(I) = 10
          if(rank.eq.0)
     &    WRITE (6,38)  NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2), SEMI1,
     &                  ECC, ECC1, HI, QPERI, R1
   38     FORMAT (' CIRCULARIZED    NAM K* AF E0 EF HI QP R* ',
     &                         2I6,2I3,1P,E10.2,0P,2F8.3,F9.1,1P,2E10.2)
*       Check time interval until Roche overflow.
          IF (KZ(34).GT.0) THEN
              CALL TRFLOW(IPAIR,DTR)
              TEV(I) = TIME + DTR
*       Update TEV but leave ROCHE call for treatment in MDOT (16/08/2006). 
*             IF (DTR.LT.STEP(I1)) THEN
*                 CALL ROCHE(IPAIR)
*                 GO TO 50
*             END IF
          END IF
      END IF
*
*       See whether a low-eccentricity synchronous state has been reached.
      RCOLL = 0.75*(RADIUS(I1) + RADIUS(I2))
      IF (ABS(SEMI1).LT.1.5*RCOLL.AND.ECC.LT.ECCM.AND.
     &    KSTAR(I).LT.10) THEN
          KSTAR(I) = 10
          if(rank.eq.0)
     &    WRITE (6,40)  ECC1, SEMI1, R(IPAIR), RCOLL, R1
   40     FORMAT (' INACTIVE PHASE    E A R RC R* ',F7.3,1P,4E10.2)
*       Check time interval until Roche overflow.
          IF (KZ(34).GT.0) THEN
              CALL TRFLOW(IPAIR,DTR)
              TEV(I) = TIME + DTR
*       Update TEV but leave ROCHE call for treatment in MDOT (16/08/2006). 
*             IF (DTR.LT.STEP(I1)) THEN
*                 CALL ROCHE(IPAIR)
*                 GO TO 50
*             END IF
          END IF
      END IF
*
*       Include warning in case of eccentricity increase (PT only).
      ECC2 = 1.0 - R(IPAIR)/SEMI1
      IF (KZ(27).EQ.1.AND.ECC2.GT.MAX(ECC,ECCM)+1.0D-04) THEN
          if(rank.eq.0)
     &    WRITE (6,42)  TTOT, IPAIR, ECC2, ECC, R(IPAIR), SEMI1
   42     FORMAT (' WARNING!    E > E0    T IP E E0 R A ',
     &                                    F10.4,I5,2F8.4,1P,2E10.2)
      END IF
*
*       Reduce radius by DR/R to delay dissipation for small energy loss.
      IF (KZ(27).EQ.1.AND.DE(1)+DE(2).LT.1.0D-07*ZMU*ABS(HI)) THEN
          J1 = I1
          IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
***       IF (TEV(J1) - TIME.GT.0.01*TEV(J1)) TEV(J1) = TIME
          DR = (0.99*4.0*RADIUS(J1) - QPERI)/QPERI
          YF = MAX(ABS(DR),0.01D0)
          RADIUS(J1) = (1.0 - MIN(YF,0.1D0))*RADIUS(J1)
          DE1 = (DE(1) + DE(2))/(ZMU*ABS(H(IPAIR)))
          if(rank.eq.0)
     &    WRITE (6,44)  TTOT, KSTAR(J1), H(IPAIR), QPERI, DR, DE1
   44     FORMAT (' REDUCED RADIUS    T K* H QP DR/R DE/E ',
     &                                F9.2,I3,F10.1,1P,3E9.1)
      END IF
*
*       Combine the two stars inelastically in case of chaos disruption.
   45 IF (KZ(27).EQ.2.AND.IDIS.GT.0) THEN
          if(rank.eq.0)
     &    WRITE (6,48) TTOT, IPAIR, LIST(1,I1), ECC, SEMI, QPERI,
     &                 RADIUS(I1), RADIUS(I2)
   48     FORMAT (' DISRUPTED CHAOS    T KS NP E A QP R* ',
     &                                 F9.2,I6,I4,F8.4,1P,4E10.2)
          CALL KSPERI(IPAIR)
          CALL XVPRED(I,0)
          KSPAIR = IPAIR
          IQCOLL = 1
          CALL CMBODY(R(IPAIR),2)
      END IF
*
   50 RETURN
*
      END
