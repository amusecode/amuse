      SUBROUTINE START4(ISUB)
*
*
*       Initialization & restart of four-body system.
*       ---------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  M,MIJ
      COMMON/CREG/  M(4),X4(3,4),XDOT4(3,4),P(12),Q(12),TIME4,ENERGY,
     &              EPSR2,XR(9),W(9),RR(6),TA(6),MIJ(6),CM(10),RMAX4,
     &              TMAX,DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/CONFIG/  R2(4,4),I1,I2,I3,I4
      COMMON/CLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3,IP(4)
      COMMON/CCOLL/  QK(12),PK(12),ICALL,ICOLL,NDISS4
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
*
*
*       Decide between new run, termination or collision (= 0, > 0, < 0).
      IF (ISUB.NE.0) THEN
          ITERM = ISUB
          ISUB = IABS(ISUB)
          GO TO 100
      END IF
*
      DO 2 K = 1,7
          CM(K) = 0.0D0
    2 CONTINUE
*
*       Transform to the local c.m. reference frame.
      DO 4 L = 1,4
          J = 2*NPAIRS + L
*       Copy four-body system from the first single particle locations.
          M(L) = BODY(J)
          CM(7) = CM(7) + M(L)
          DO 3 K = 1,3
              X4(K,L) = X(K,J)
              XDOT4(K,L) = XDOT(K,J)
              CM(K) = CM(K) + M(L)*X4(K,L)
              CM(K+3) = CM(K+3) + M(L)*XDOT4(K,L)
    3     CONTINUE
    4 CONTINUE
*
*       Set c.m. coordinates & velocities of four-body system.
      DO 5 K = 1,6
          CM(K) = CM(K)/CM(7)
    5 CONTINUE
*
*       Specify initial conditions for chain regularization.
      DO 8 L = 1,4
          DO 7 K = 1,3
              X4(K,L) = X4(K,L) - CM(K)
              XDOT4(K,L) = XDOT4(K,L) - CM(K+3)
    7     CONTINUE
    8 CONTINUE
*
*       Calculate internal energy and include in total subsystem energy.
      CALL NEWSYS(X4,XDOT4,M,4,ENERGY,GAM)
      ESUB = ESUB + ENERGY
*
*       Save global indices & particle attributes of subsystem.
      DO 10 L = 1,4
          J = 2*NPAIRS + L
          JLIST(L) = J
          NAME4(L) = NAME(J)
          SIZE(L) = RADIUS(J)
          IP(L) = KSTAR(J)
   10 CONTINUE
*
*       Make perturber list for potential energy correction.
      ICM = JLIST(1)
      GO TO 200
*
*       Obtain potential enegy of resolved subsystem & perturbers.
   20 CALL NBPOT(4,NP,POT1)
*
*       Define subsystem indicator (ISYS = 1, 2, 3 for triple, quad, chain).
      ISYS(NSUB+1) = 2
*
*       Form ghosts and initialize c.m. motion in ICOMP (= JLIST(1)).
      CALL SUBSYS(4,CM)
*
*       Include interaction of subsystem c.m. & perturbers for net effect.
      CALL NBPOT(1,NP,POT2)
*
*       Form square of c.m. velocity correction due to differential force.
      VI2 = X0DOT(1,ICOMP)**2 + X0DOT(2,ICOMP)**2 + X0DOT(3,ICOMP)**2
      CORR = 1.0 + 2.0*(POT2 - POT1)/(CM(7)*VI2)
      IF (CORR.LE.0.0D0) CORR = 0.0
*
*       Modify c.m. velocity by net tidal energy correction.
      DO 30 K = 1,3
          X0DOT(K,ICOMP) = SQRT(CORR)*X0DOT(K,ICOMP)
   30 CONTINUE
*
*       Remove ghosts from perturber neighbour lists.
      CALL NBREM(ICM,4,NP)
*
*       Set maximum integration interval equal to c.m. step.
      TMAX = STEP(ICOMP)
*
*       Copy total energy and output & capture options for routine QUAD.
      CM(8) = BE(3)
      KZ15 = KZ(15)
      KZ27 = KZ(27)
*
*       Assign new subsystem index and begin four-body regularization.
      ISUB = NSUB
      NQUAD = NQUAD + 1
      GO TO 180
*
*       Prepare KS regularization and direct integration of two bodies.
  100 JLIST(5) = NAME4(I1)
      JLIST(6) = NAME4(I2)
      JLIST(7) = NAME4(I3)
      JLIST(8) = NAME4(I4)
*
*       Indentify current global index by searching all single particles.
      DO 102 J = IFIRST,N
          DO 101 L = 1,4
              IF (NAME(J).EQ.JLIST(L+4)) THEN
                  JLIST(L) = J
              END IF
  101     CONTINUE
  102 CONTINUE
*
*       Ensure ICOMP < JCOMP for KS regularization.
      ICOMP = MIN(JLIST(1),JLIST(2))
      JCOMP = MAX(JLIST(1),JLIST(2))
*
*       Identify global index of c.m. body.
      ICM = 0
      DO 104 L = 1,4
          J = JLIST(L)
          IF (BODY(J).GT.0.0D0) ICM = J
  104 CONTINUE
*
*       Quantize the elapsed interval since last step.
      TIME2 = T0S(ISUB) + TIME4 - TPREV
      DT8 = (TBLOCK - TPREV)/8.0D0
*
*       Adopt the nearest truncated step (at most 8 subdivisions).
      DT2 = TIME2
      IF (TIME2.GT.0.0D0) THEN
          CALL STEPK(DT2,DTN2)
          DTN = NINT(DTN2/DT8)*DT8
      ELSE
*       Choose negative step if pericentre time < TPREV (cf. iteration).
          DT2 = -DT2
          CALL STEPK(DT2,DTN2)
          DTN = -NINT(DTN2/DT8)*DT8
      END IF
*
*       Update time for new polynomial initializations (also for CMBODY).
      TIME = TPREV + DTN
      TIME = MIN(TBLOCK,TIME)
*
*       Predict current coordinates & velocities before termination.
      CALL XVPRED(ICM,0)
*
*       Coopy c.m. coordinates & velocities.
      DO 105 K = 1,3
          CM(K) = X(K,ICM)
          CM(K+3) = XDOT(K,ICM)
  105 CONTINUE
*
*       Re-determine the perturber list.
      GO TO 200
*
*       Obtain potential energy of the c.m. subsystem & JPERT(NP).
  110 I = JLIST(1)
      JLIST(1) = ICM
      CALL NBPOT(1,NP,POT1)
*
*       Set configuration pointers for KS candidates & distant bodies.
      JLIST(5) = I1
      JLIST(6) = I2
      JLIST(7) = I3
      JLIST(8) = I4
*
*       Place new coordinates in the original locations.
      JLIST(1) = I
      DO 120 L = 1,4
          J = JLIST(L)
*       Compare global name & subsystem name to restore the mass.
          DO 112 K = 1,4
              IF (NAME(J).EQ.NAMES(K,ISUB)) THEN
                  BODY(J) = BODYS(K,ISUB)
              END IF
  112     CONTINUE
          LL = JLIST(L+4)
          DO 115 K = 1,3
              X(K,J) = X4(K,LL) + CM(K)
  115     CONTINUE
  120 CONTINUE
*
*       Obtain potential energy of subsystem & perturbers at the end.
      CALL NBPOT(4,NP,POT2)
*
*       Form square of c.m. velocity correction due to differential force.
      VI2 = CM(4)**2 + CM(5)**2 + CM(6)**2
      CORR = 1.0 + 2.0*(POT2 - POT1)/(CM(7)*VI2)
      IF (CORR.LE.0.0D0) CORR = 0.0
*
*       Modify c.m. velocity by net tidal energy correction.
      DO 122 K = 1,3
          CM(K+3) = SQRT(CORR)*CM(K+3)
  122 CONTINUE
*
*       Transform to global velocities using corrected c.m. velocity.
      DO 130 L = 1,4
          J = JLIST(L)
          LL = JLIST(L+4)
          DO 125 K = 1,3
              XDOT(K,J) = XDOT4(K,LL) + CM(K+3)
              X0DOT(K,J) = XDOT(K,J)
  125     CONTINUE
  130 CONTINUE
*
*       Predict coordinates & velocities of perturbers to order FDOT.
      DO 140 L = 1,NP
          J = JPERT(L)
          CALL XVPRED(J,0)
  140 CONTINUE
*
*       Update subsystem COMMON variables unless last or only case.
      IF (ISUB.LT.NSUB) THEN
          DO 150 L = ISUB,NSUB
              DO 145 K = 1,6
                  BODYS(K,L) = BODYS(K,L+1)
                  NAMES(K,L) = NAMES(K,L+1)
  145         CONTINUE
              T0S(L) = T0S(L+1)
              TS(L) = TS(L+1)
              STEPS(L) = STEPS(L+1)
              RMAXS(L) = RMAXS(L+1)
              ISYS(L) = ISYS(L+1)
  150     CONTINUE
      END IF
*
*       Reduce subsystem counter and subtract internal binding energy.
      NSUB = NSUB - 1
      ESUB = ESUB - ENERGY - ECOLL3
*
*       Select neighbours for single bodies and new c.m. (JLIST is *2).
      RS0 = RS(ICM)
      J3 = JLIST(3)
      J4 = JLIST(4)
      CALL NBLIST(J3,RS0)
      CALL NBLIST(J4,RS0)
      CALL NBLIST(ICOMP,RS0)
*
*       Add any other ghosts to perturber list for replacement of #ICM.
      DO 160 JSUB = 1,NSUB
          DO 155 L = 1,4
          IF (NAMES(L,JSUB).EQ.0) GO TO 155
              DO 154 J = 1,N
              IF (NAME(J).EQ.NAMES(L,JSUB).AND.BODY(J).LE.0.0D0) THEN
                  NP = NP + 1
                  JPERT(NP) = J
                  GO TO 155
              END IF
  154         CONTINUE
  155     CONTINUE
  160 CONTINUE
*
*       Replace ICM in perturber neighbour lists by all subsystem members.
      CALL NBREST(ICM,4,NP)
*
*       Check for stellar collision (only needs coordinates & velocities).
      IF (ITERM.LT.0) THEN
          JLIST(1) = ICOMP
          JLIST(2) = JCOMP
*
*       See whether relabelling is required (indices I1 - I4 still local).
          IF (R2(I1,I4).LT.R2(I1,I3).OR.R2(I3,I4).LT.R2(I1,I3)) THEN
              IF (R2(I1,I4).LT.R2(I3,I4)) THEN
*       Switch body #I3 & I4 to give new dominant pair I1 & I3.
                  I = JLIST(4)
                  JLIST(4) = JLIST(3)
                  JLIST(3) = I
              ELSE
*       Set JLIST(5) < 0 to denote that body #I3 & I4 will be new KS pair.
                  JLIST(5) = -1
              END IF
          END IF
          GO TO 170
      END IF
*
*       Exclude the dominant interaction for c.m. approximation (large FDOT).
      IF (MIN(R2(I1,I3),R2(I2,I4)).GT.CMSEP2*R2(I1,I2)) THEN
          JLIST(1) = JLIST(I3)
          JLIST(2) = JLIST(I4)
          NNB = 2
      ELSE
          NNB = 4
      END IF
*
*       Set dominant F & FDOT on body #ICOMP & JCOMP for #I3 & I4 in FPOLY2.
      CALL FCLOSE(ICOMP,NNB)
      CALL FCLOSE(JCOMP,NNB)
*
*       Specify global indices of least dominant bodies.
      I3 = JLIST(3)
      I4 = JLIST(4)
*
*       Initialize force polynomials & time-steps for body #I3 & #I4.
      CALL FPOLY1(I3,I3,0)
      CALL FPOLY1(I4,I4,0)
*
      CALL FPOLY2(I3,I3,0)
      CALL FPOLY2(I4,I4,0)
*
*       Perform KS regularization of dominant components (ICOMP < JCOMP).
      CALL KSREG
*
*       Check minimum two-body distance.
      DMIN4 = MIN(DMIN4,RCOLL)
*
*       Update net binary energy change.
      BBCOLL = BBCOLL + CM(9)
*
*       Update number of DIFSY calls, tidal dissipations & collision energy.
  170 NSTEPQ = NSTEPQ + NSTEP4
      NDISS = NDISS + NDISS4
      ECOLL = ECOLL + ECOLL3
      E(10) = E(10) + ECOLL3
*
*       Check for subsystem at last COMMON dump (no restart with NSUB > 0).
      IF (NSUB.EQ.0.AND.KZ(2).GE.1) THEN
          IF (TIME - TDUMP.LT.TIME4) THEN
              TDUMP = TIME
              CALL MYDUMP(1,2)
          END IF
      END IF
*
*       Set phase indicator = -1 to ensure new time-step list in INTGRT.
  180 IPHASE = -1
*
      RETURN
*
*       Form the current perturber list.
  200 RP2 = RS(ICM)**2
      NP = 0
      NNB2 = LIST(1,ICM) + 1
*
*       Loop over all single particles & c.m. but skip subsystem members.
      DO 210 L = 2,NNB2
          J = LIST(L,ICM)
          RIJ2 = (X(1,J) - CM(1))**2 + (X(2,J) - CM(2))**2 +
     &                                 (X(3,J) - CM(3))**2
          IF (RIJ2.LT.RP2) THEN
              DO 205 K = 1,4
                  IF (J.EQ.JLIST(K)) GO TO 210
  205         CONTINUE
              NP = NP + 1
              JPERT(NP) = J
          END IF
  210 CONTINUE
*
      IF (ISUB.EQ.0) GO TO 20
      GO TO 110
*
      END
