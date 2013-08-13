      SUBROUTINE CHECKL(I,NNB,XI,XIDOT,RS2,FIRR,FREG,FD,FDR)
*
*
*       Neighbour list modification.
*       ----------------------------
*
      INCLUDE 'common6.h'
      REAL*8  XI(3),XIDOT(3),DX(3),DV(3),FIRR(3),FREG(3),FD(3),FDR(3)
*
*
      ICM = 1
*       Only consider high-velocity particles if KZ(37) = 1.
      IF (KZ(37).EQ.1) GO TO 350
*       Omit special treatments of c.m. particles because of force errors.
      IF (I.GT.N) GO TO 350
*
*       See whether any neighbour has encounter with body outside sphere.
      NNB1 = NNB + 1
      LJ = 0
*       First perform fast test since actual cases are rare.
      DO 306 L = 2,NNB1
          J = ILIST(L)
          IF (STEP(J).LT.SMIN) LJ = L
  306 CONTINUE
      IF (LJ.EQ.0) GO TO 320
*
      JCOMP = ILIST(LJ)
      K = 1
  308 IF (K.GT.LIST(1,JCOMP)) GO TO 320
      K = K + 1
      J = LIST(K,JCOMP)
      IF (STEP(J).GT.SMIN) GO TO 308
*       Skip single regularized particles or body #I itself.
      IF (J.LT.IFIRST.OR.J.EQ.I) GO TO 308
      A1 = X(1,J) - X(1,JCOMP)
      A2 = X(2,J) - X(2,JCOMP)
      A3 = X(3,J) - X(3,JCOMP)
      RIJ2 = A1**2 + A2**2 + A3**2
*       Only accept body #J as neighbour if distance to JCOMP is < 2*RMIN.
      IF (RIJ2.GT.RMIN22) GO TO 308
*
      IF (J.LE.N) GO TO 309
*       Also accept pairs satisfying the c.m. force approximation.
      IF (CMSEP2*R(J-N)**2.GT.RS2) GO TO 308
*
*       Now see whether body #J has already been included as neighbour.
  309 DO 310 L = 1,NNB
          IF (J.EQ.ILIST(L+1)) GO TO 308
  310 CONTINUE
*       Include body #J in sequential location of neighbour list.
      L = NNB + 1
  312 IF (J.GT.ILIST(L)) GO TO 314
      ILIST(L+1) = ILIST(L)
*       Move other members down by one until appropriate location is free.
      L = L - 1
      IF (L.GT.1) GO TO 312
*
  314 ILIST(L+1) = J
      NNB = NNB + 1
*       Only one close encounter neighbour is allowed without test of NNB.
      NLSMIN = NLSMIN + 1
*
*       Finally correct irregular and regular force components.
      DR2 = 0.0
      DRDV = 0.0
      DO 315 K = 1,3
          DX(K) = X(K,J) - XI(K)
          DV(K) = XDOT(K,J) - XIDOT(K)
          DR2 = DR2 + DX(K)**2
          DRDV = DRDV + DX(K)*DV(K)
  315 CONTINUE
*
      DR2I = 1.0/DR2
      DR3I = BODY(J)*DR2I*SQRT(DR2I)
      DRDV = 3.0*DRDV*DR2I
*
      DO 316 K = 1,3
          FIRR(K) = FIRR(K) + DX(K)*DR3I
          FREG(K) = FREG(K) - DX(K)*DR3I
          FD(K) = FD(K) + (DV(K) - DX(K)*DRDV)*DR3I
          FDR(K) = FDR(K) - (DV(K) - DX(K)*DRDV)*DR3I
  316 CONTINUE
*
*       Include the other component of two recently disrupted pairs.
  320 IF (NNB.GE.NNBMAX) GO TO 330
*
      L = 0
  321 L = L + 2
*       Advance list index by two for every identified pair.
      IF (ILIST(L).GT.IFIRST + 3.OR.L.GT.NNB + 1) GO TO 330
*
*       Set appropriate pair index for either component.
      JL = ILIST(L)
      JPAIR = KVEC(JL)
      IL1 = ILIST(L+1)
      IPAIR = KVEC(IL1)
*       Check whether two consecutive list members belong to same pair.
      IF (JPAIR.EQ.IPAIR) GO TO 321
*       The case of index L referring to the last neighbour is permitted.
      J = 2*JPAIR
      IF (J.EQ.ILIST(L)) J = J - 1
*       Index of the missing component, subject to neighbour test.
      IF (J.EQ.I) GO TO 323
      JCOMP = ILIST(L)
      A1 = X(1,JCOMP) - X(1,J)
      A2 = X(2,JCOMP) - X(2,J)
      A3 = X(3,JCOMP) - X(3,J)
      RIJ2 = A1*A1 + A2*A2 + A3*A3
*       Only accept #J as neighbour if distance to JCOMP is < 2*RMIN.
      IF (RIJ2.LT.RMIN22) GO TO 324
  323 L = L - 1
*       Increase search index by one only after unsuccessful test.
      GO TO 321
*
*       Correct irregular & regular force components.
  324 DR2 = 0.0
      DRDV = 0.0
      DO 325 K = 1,3
          DX(K) = X(K,J) - XI(K)
          DV(K) = XDOT(K,J) - XIDOT(K)
          DR2 = DR2 + DX(K)**2
          DRDV = DRDV + DX(K)*DV(K)
  325 CONTINUE
*
      DR2I = 1.0/DR2
      DR3I = BODY(J)*DR2I*SQRT(DR2I)
      DRDV = 3.0*DRDV*DR2I
*
      DO 326 K = 1,3
          FIRR(K) = FIRR(K) + DX(K)*DR3I
          FREG(K) = FREG(K) - DX(K)*DR3I
          FD(K) = FD(K) + (DV(K) - DX(K)*DRDV)*DR3I
          FDR(K) = FDR(K) - (DV(K) - DX(K)*DRDV)*DR3I
  326 CONTINUE
*
      LJ = NNB + 1
  327 IF (J.GT.ILIST(LJ)) GO TO 328
*       Move other members down by one until relevant location is free.
      ILIST(LJ+1) = ILIST(LJ)
      LJ = LJ - 1
      IF (LJ.GT.1) GO TO 327
*
  328 ILIST(LJ+1) = J
      NNB = NNB + 1
      NBDIS = NBDIS + 1
*       Note that next list index is increased by two after including #J.
      IF (NNB.LT.NNBMAX) GO TO 321
*
*       This part includes the other component of an exchanged pair.
  330 IF (LISTR(1).EQ.0.OR.NNB.GE.NNBMAX) GO TO 350
*
*       Do a fast skip if only one disrupted pair in original location.
      IF (LISTR(1).EQ.2.AND.LISTR(2).LE.IFIRST + 3) GO TO 350
      NNB1 = LISTR(1)
      LJ = 1 + 0.2*FLOAT(NNB)
      L = 1
  332 ICASE = 0
      LG = 0
      KTIME = 0
  334 KTIME = KTIME + 1
  335 IF (L.GT.NNB1) GO TO 350
*
      L = L + 1
      JBODY = LISTR(L)
      IF (JBODY.LE.IFIRST + 3) GO TO 335
*       The two last disrupted pairs have already been considered.
  336 LG = LG + LJ
*       First use large increments to speed up the search.
      IF (LG.GT.NNB + 1) LG = NNB + 1
      IF (ILIST(LG).LT.JBODY.AND.LG.LT.NNB + 1) GO TO 336
*
      LG = LG + 1
  338 LG = LG - 1
      IF (ILIST(LG).GT.JBODY.AND.LG.GT.2) GO TO 338
      IF (ILIST(LG).EQ.JBODY) ICASE = ICASE + KTIME
      IF (KTIME.EQ.1) GO TO 334
*
*       See whether any or both of the components have been identified.
      IF (ICASE.EQ.0.OR.ICASE.EQ.3) GO TO 332
*
      J = LISTR(L+1-ICASE)
*       Index of missing component to be included subject to J = I test.
      JCOMP = LISTR(L+ICASE-2)
*       Arguments for J & JCOMP are L or L - 1 depending on ICASE.
      A1 = X(1,JCOMP) - X(1,J)
      A2 = X(2,JCOMP) - X(2,J)
      A3 = X(3,JCOMP) - X(3,J)
      RIJ2 = A1*A1 + A2*A2 + A3*A3
*       Accept body #J only if distance to JCOMP is < 2*RMIN (skip J = I).
      IF (RIJ2.GT.RMIN22.OR.J.EQ.I) GO TO 332
*
*       Carry out force modifications due to addition of neighbour.
  342 DR2 = 0.0
      DRDV = 0.0
      DO 343 K = 1,3
          DX(K) = X(K,J) - XI(K)
          DV(K) = XDOT(K,J) - XIDOT(K)
          DR2 = DR2 + DX(K)**2
          DRDV = DRDV + DX(K)*DV(K)
  343 CONTINUE
*
      DR2I = 1.0/DR2
      DR3I = BODY(J)*DR2I*SQRT(DR2I)
      DRDV = 3.0*DRDV*DR2I
*
      DO 344 K = 1,3
          FIRR(K) = FIRR(K) + DX(K)*DR3I
          FREG(K) = FREG(K) - DX(K)*DR3I
          FD(K) = FD(K) + (DV(K) - DX(K)*DRDV)*DR3I
          FDR(K) = FDR(K) - (DV(K) - DX(K)*DRDV)*DR3I
  344 CONTINUE
*
*       Include body #J in neighbour list and increase NNB.
  345 K = NNB + 1
  346 IF (J.GT.ILIST(K)) GO TO 348
*       Move other members down by one until relevant location is free.
      ILIST(K+1) = ILIST(K)
      K = K - 1
      IF (K.GT.1) GO TO 346
*
  348 ILIST(K+1) = J
      NNB = NNB + 1
      NBDIS2 = NBDIS2 + 1
*
*       Check list of high velocity intruders for early retention.
  350 IF (LISTV(1).EQ.0.OR.ICM.EQ.-1) GO TO 355
      L = 2
  351 J = LISTV(L)
      IF (J.EQ.I) GO TO 354
      A1 = X(1,J) - XI(1)
      A2 = X(2,J) - XI(2)
      A3 = X(3,J) - XI(3)
      RIJ2 = A1**2 + A2**2 + A3**2
*     IF (RIJ2.GT.4.0*RS2.OR.RIJ2.LT.RS2) GO TO 354
*       Simplify to include GPU (standard code needs a few extra searches).
      IF (RIJ2.GT.4.0*RS2) GO TO 354
      A4 = XDOT(1,J) - XDOT(1,I)
      A5 = XDOT(2,J) - XDOT(2,I)
      A6 = XDOT(3,J) - XDOT(3,I)
      A7 = A1*A4 + A2*A5 + A3*A6
      IF (A7.GT.0.0) GO TO 354
      P2 = RIJ2 - A7**2/(A4**2 + A5**2 + A6**2)
*       Accept if impact parameter < RS.
      IF (P2.GT.RS2) GO TO 354
*       See if body #J has been included by other procedures.
      DO 353 K = 1,NNB
          IF (J.EQ.ILIST(K+1)) GO TO 354
  353 CONTINUE
      ICM = -1
*       Redundant indicator (normally > 0) used for joint procedure.
      NBDIS2 = NBDIS2 - 1
      NBFAST = NBFAST + 1
*
*       Distinguish betweeen single particle treatment and perturbed case.
      IF (NNB.LE.NNBMAX) THEN
          IF (I.LE.N) GO TO 342
          RIJ2 = (XI(1) - X(1,J))**2 + (XI(2) - X(2,J))**2 +
     &                                 (XI(3) - X(3,J))**2
*       Adopt c.m. approximation if body #J is also single.
          IF (RIJ2.GT.CMSEP2*R(I-N)**2.AND.J.LE.N) GO TO 342
      ELSE
          GO TO 355
      END IF 
*
*       Set KS component indices (also body #J if inside c.m. approximation).
      I1 = 2*(I - N) - 1
      I2 = I1 + 1
      J1 = J
      J2 = 0
      IF (J.GT.N) THEN
          IF (RIJ2.LT.CMSEP2*R(J-N)**2) THEN
              J1 = 2*(J - N) - 1
              J2 = J1 + 1
          END IF
      END IF
*
*       Begin evaluation of all relevant interactions (at most 4 terms).
      K = J1
  360 L = I1
  362 A1 = X(1,K) - X(1,L)
      A2 = X(2,K) - X(2,L)
      A3 = X(3,K) - X(3,L)
      RIJ2 = A1*A1 + A2*A2 + A3*A3
      DV(1) = XDOT(1,K) - XDOT(1,L)
      DV(2) = XDOT(2,K) - XDOT(2,L)
      DV(3) = XDOT(3,K) - XDOT(3,L)
*
*       Employ the appropriate mass-weighted expression.
      DR2I = 1.0/RIJ2
      DR3I = BODY(K)*BODY(L)*DR2I*SQRT(DR2I)/BODY(I)
      DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
*       Add irregular F & FDOT and subtract regular terms.
      FIRR(1) = FIRR(1) + A1*DR3I
      FIRR(2) = FIRR(2) + A2*DR3I
      FIRR(3) = FIRR(3) + A3*DR3I
      FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
      FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
      FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
      FREG(1) = FREG(1) - A1*DR3I
      FREG(2) = FREG(2) - A2*DR3I
      FREG(3) = FREG(3) - A3*DR3I
      FDR(1) = FDR(1) - (DV(1) - A1*DRDV)*DR3I
      FDR(2) = FDR(2) - (DV(2) - A2*DRDV)*DR3I
      FDR(3) = FDR(3) - (DV(3) - A3*DRDV)*DR3I
*
*       Consider each interaction in turn (body #J may be resolved).
      L = L + 1
      IF (L.EQ.I2) GO TO 362
      K = K + 1
      IF (K.EQ.J2) GO TO 360
*
*       Include body #J in the neighbour list.
      GO TO 345
*
  354 L = L + 1
      IF (L.LE.LISTV(1) + 1) GO TO 351
*
  355 RETURN
*
      END
