      SUBROUTINE KSINIT
*
*
*       Initialization of KS regularization.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  Q(3),RDOT(3),UI(4),VI(4),A1(3,4)
      SAVE  IPREV,EBPREV
      DATA  IPREV,EBPREV /0,-1.0D-06/
*
*
*       Set new global indices of the components and current pair index.
      ICOMP = 2*NPAIRS - 1
      JCOMP = ICOMP + 1
      IPAIR = NPAIRS
*
*       Add body #N in case the only neighbour was removed in KSREG.
      IF (LIST(1,ICOMP).EQ.0) THEN
          LIST(2,ICOMP) = N
          LIST(2,JCOMP) = N
          LIST(1,ICOMP) = 1
          LIST(1,JCOMP) = 1
      END IF
*
*       Specify mass, neighbour radius, name, radius & type for new c.m.
      BODY(NTOT) = BODY(ICOMP) + BODY(JCOMP)
      RS(NTOT) = RS(ICOMP)
      NAME(NTOT) = NZERO + NAME(ICOMP)
      RADIUS(NTOT) = 0.0
      TEV(NTOT) = 1.0E+10
      TEV0(NTOT) = 1.0E+10
      BODY0(NTOT) = BODY(NTOT)
      EPOCH(NTOT) = TIME*TSTAR
      KSTAR(NTOT) = 0
      IMOD = 1
*
*       Define c.m. coordinates & velocities and set XDOT for components.
      DO 10 K = 1,3
          X(K,NTOT) = (BODY(ICOMP)*X(K,ICOMP) + BODY(JCOMP)*X(K,JCOMP))/
     &                                                        BODY(NTOT)
          X0DOT(K,NTOT) = (BODY(ICOMP)*X0DOT(K,ICOMP) + BODY(JCOMP)*
     &                                        X0DOT(K,JCOMP))/BODY(NTOT)
          XDOT(K,NTOT) = X0DOT(K,NTOT)
          XDOT(K,ICOMP) = X0DOT(K,ICOMP)
          XDOT(K,JCOMP) = X0DOT(K,JCOMP)
   10 CONTINUE
*
*       Obtain force polynomial for c.m. with components ICOMP & JCOMP.
      NNB = LIST(1,ICOMP)
*
*       Predict current coordinates & velocities for the neighbours.
      CALL XVPRED(ICOMP,NNB)
*
*       Obtain new polynomials & steps (first F & FDOT, then F2DOT & F3DOT).
      CALL FPOLY1(ICOMP,JCOMP,1)
      CALL FPOLY2(NTOT,NTOT,1)
*
*       Skip KS initialization at merger termination (H, U & UDOT in RESET).
      IF (IPHASE.EQ.7) THEN
          EB = 2.0*EBH
          GO TO 50
      END IF
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
*       Form regularized velocity and set initial KS coordinates.
      TDOT2(IPAIR) = 0.0
      DO 30 K = 1,4
          UDOT(K,IPAIR) = 0.50D0*(A1(1,K)*RDOT(1) + A1(2,K)*RDOT(2) +
     &                                              A1(3,K)*RDOT(3))
*       Note that A1(J,K) is the transpose of A1(K,J).
          U(K,IPAIR) = UI(K)
          U0(K,IPAIR) = U(K,IPAIR)
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0D0*UI(K)*UDOT(K,IPAIR)
   30 CONTINUE
*
*       Evaluate initial binding energy per unit mass.
      H(IPAIR) = (2.0D0*(UDOT(1,IPAIR)**2 + UDOT(2,IPAIR)**2 +
     &                   UDOT(3,IPAIR)**2 + UDOT(4,IPAIR)**2) -
     &                                              BODY(NTOT))/R(IPAIR)
      EB = H(IPAIR)*BODY(ICOMP)*BODY(JCOMP)/BODY(NTOT)
*
*       Form perturber list.
   50 CALL KSLIST(IPAIR)
*
*       Transform any unperturbed hard binary to apocentre and set time-step.
      SEMI = -0.5*BODY(NTOT)/H(IPAIR)
      IF (LIST(1,ICOMP).EQ.0.AND.EB.LT.EBH) THEN
          TK = TWOPI*SEMI*SQRT(SEMI/BODY(NTOT))
*       Note TIME is not commensurate after KSPERI (cf. CHTERM & STEPS).
          IF (IPHASE.NE.7.AND.IPHASE.NE.8) THEN
              DO 55 K = 1,4
                  VI(K) = UDOT(K,IPAIR)
   55         CONTINUE
*       Determine pericentre time (TP < 0 if TDOT2 < 0) and add TK/2.
              CALL TPERI(SEMI,UI,VI,BODY(NTOT),TP)
*       Note: apocentre to apocentre gives almost zero step.
              STEP(ICOMP) = 0.5*MIN(TK,STEP(NTOT)) - TP
*       Transform KS variables to peri and by pi/2 to apocentre (skip apo).
              IF (ABS(TDOT2(IPAIR)).GT.1.0E-12.OR.R(IPAIR).LT.SEMI) THEN
                  TIME0 = TIME
                  CALL KSPERI(IPAIR)
                  CALL KSAPO(IPAIR)
*       Reset TIME to quantized value (small > 0 or < 0 possible initially).
                  TIME = TIME0
              ELSE IF (TDOT2(IPAIR).GT.0.0) THEN
                  TDOT2(IPAIR) = -1.0E-20
              END IF
          END IF
      END IF
*
*       Estimate an appropriate KS slow-down index for G < GMIN.
      IF (LIST(1,ICOMP).EQ.0.AND.SEMI.GT.0.0) THEN
          TK = TWOPI*SEMI*SQRT(SEMI/BODY(NTOT))
          IF (KZ(26).GT.0.AND.STEP(NTOT).GT.TK) THEN
              IMOD = 1 + LOG(STEP(NTOT)/TK)/0.69
              IMOD = MIN(IMOD,5)
          END IF
      END IF
*
*       Specify zero membership and large step for second component).
      LIST(1,JCOMP) = 0
*       Set large step for second component to avoid detection.
      STEP(JCOMP) = 1.0E+06
      STEPR(JCOMP) = 1.0E+06
*
*       Obtain polynomials for perturbed KS motion (standard case & merger).
      CALL KSPOLY(IPAIR,IMOD)
*
*       Obtain apocentre distance.
      SEMI = -0.5*BODY(NTOT)/H(IPAIR)
      ECC2 = (1.0-R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(NTOT)*SEMI)
      RAP = SEMI*(1.0 + SQRT(ECC2))
*
*       Include suggestion for monitoring hyperbolic encounters (suppressed).
*     IF (SEMI.LT.0.0) THEN
*         PMIN = SEMI*(1.0 - SQRT(ECC2))
*         if(rank.eq.0)
*    &    WRITE (31,56)  TIME+TOFF, NAME(ICOMP), NAME(JCOMP), PMIN
*  56     FORMAT (' HYPERBOLIC    T NAM PMIN ',F8.2,2I6,1P,E10.2)
*     END IF
*
*       Set 2*SEMI as termination scale for hard binary if 2*SEMI < RS/20.
      IF (EB.LT.EBH.AND.RAP.LT.0.05*RS(NTOT)) THEN 
          R0(IPAIR) = MAX(RMIN,-BODY(NTOT)/H(IPAIR))
      ELSE
          R0(IPAIR) = R(IPAIR)
      END IF
*
*       Increase regularization counters (#9 & NKSHYP for hyperbolic orbits).
      NKSREG = NKSREG + 1
      IF (H(IPAIR).GT.0.0) NKSHYP = NKSHYP + 1
*
*       Check optional output for new KS.
      IF (KZ(10).GT.0) THEN
          RI = SQRT((X(1,NTOT) - RDENS(1))**2 +
     &              (X(2,NTOT) - RDENS(2))**2 +
     &              (X(3,NTOT) - RDENS(3))**2)
          if(rank.eq.0)
     &    WRITE (6,60)  TIME+TOFF, NAME(ICOMP), NAME(JCOMP),DTAU(IPAIR),
     &             R(IPAIR), RI, H(IPAIR), IPAIR, NAME(N+IPAIR),
     &             GAMMA(IPAIR),STEP(NTOT), LIST(1,ICOMP), LIST(1,NTOT)
   60     FORMAT (/,' NEW KSREG    TIME =',F8.2,2I8,F12.3,1PE10.1,
     &                                0PF7.2,F9.2,2I8,F8.3,1PE10.1,2I5)
      END IF
*
*       Include diagnostics for NS or BH hard binary formation.
      IF (MAX(KSTAR(ICOMP),KSTAR(JCOMP)).GE.13.AND.EB.LT.EBH.AND.
     &    IPHASE.NE.7) THEN
*       Limit the diagnostics to significant change or new case.
          ISUM = KSTAR(ICOMP) + KSTAR(JCOMP) + KSTAR(NTOT)
          DEB = ABS((EB - EBPREV)/EBPREV)
          IF (ISUM.NE.IPREV.OR.DEB.GT.0.1) THEN
              IPREV = ISUM
              EBPREV = EB
              PD = DAYS*SEMI*SQRT(SEMI/BODY(NTOT))
              if(rank.eq.0)
     &        WRITE (6,65)  TIME+TOFF, NAME(ICOMP), NAME(JCOMP),
     &          NAME(N+IPAIR), KSTAR(ICOMP), KSTAR(JCOMP), KSTAR(NTOT),
     &                      SQRT(ECC2), PD, EB
   65         FORMAT (' NS/BH BINARY    T NM K* E P EB ',1P,
     &                                  E12.5,3I8,3I4,3E12.5)
          END IF
      END IF
*
*       Modify the termination criterion according to value of NPAIRS.
      IF (NPAIRS.GT.KMAX - 3) GMAX = 0.8*GMAX
      IF (NPAIRS.LT.KMAX - 5.AND.GMAX.LT.0.001) GMAX = 1.2*GMAX
      IF (NPAIRS.EQ.KMAX.and.rank.eq.0) WRITE (6,70)  NPAIRS, TIME+TOFF
   70 FORMAT (5X,'WARNING!   MAXIMUM KS PAIRS   NPAIRS TIME',I5,F9.2)
*
*       Initialize prediction variables to avoid skipping KS components.
      DO 75 KCOMP = 1,2
          JDUM = 2*NPAIRS - 2 + KCOMP
          DO 72 K = 1,3
              X0(K,JDUM) = X(K,JDUM)
              X0DOT(K,JDUM) = 0.0D0
              F(K,JDUM) = 0.0D0
              FDOT(K,JDUM) = 0.0D0
              D2(K,JDUM) = 0.0D0
              D3(K,JDUM) = 0.0D0
              D2R(K,JDUM) = 0.0D0
              D3R(K,JDUM) = 0.0D0
   72     CONTINUE
   75 CONTINUE
*
*       See whether either component has been regularized recently.
      NNB = LISTD(1) + 1
      K = 0
*       Check case of initial binary and loop over disrupted pairs.
      IF (IABS(NAME(ICOMP) - NAME(JCOMP)).EQ.1) THEN
          IF (NAME(ICOMP).LE.2*NBIN0) K = -1
      END IF
      DO 80 L = 2,NNB
          IF (NAME(ICOMP).EQ.LISTD(L).OR.NAME(JCOMP).EQ.LISTD(L)) K = -1
   80 CONTINUE
      IF (H(IPAIR).GT.0.0) K = 0
*
*       Treat mergers as new binaries and ensure chain/hard KS as standard.
      IF (IPHASE.EQ.6) THEN
          K = 0
      ELSE IF (K.EQ.0) THEN
          IF (IPHASE.EQ.8.OR.EB.LT.EBH) K = -1
      END IF
*       Set flag to distinguish between existing and new binaries (K = -1/0).
      LIST(2,JCOMP) = K
      KSLOW(IPAIR) = 1
*
*       Check diagnostic output of new hard binary.
      IF (KZ(8).GT.0.AND.K.EQ.0) THEN
          IF (EB.GT.EBH) GO TO 100
          SEMI = -0.5*BODY(NTOT)/H(IPAIR)
          RI = SQRT((X(1,NTOT) - RDENS(1))**2 +
     &              (X(2,NTOT) - RDENS(2))**2 +
     &              (X(3,NTOT) - RDENS(3))**2)
          if(rank.eq.0)
     &    WRITE (8,90)  TIME+TOFF, NAME(ICOMP), NAME(JCOMP), K,
     &                  BODY(ICOMP),BODY(JCOMP), EB, SEMI, R(IPAIR),
     &                  GAMMA(IPAIR), RI
   90     FORMAT (' NEW BINARY   T =',F7.1,'  NAME = ',2I6,I3,
     &                        '  M =',1P,2E9.1,'  EB =',E9.1,
     &                        '  A =',E9.1,'  R =',E9.1,'  G =',E9.1,
     &                        '  RI =',E9.1)
          CALL FLUSH(8)
      END IF
*
  100 RETURN
*
      END
