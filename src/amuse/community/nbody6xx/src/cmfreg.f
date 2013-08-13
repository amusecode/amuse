      SUBROUTINE CMFREG(I,RS2,RCRIT2,VRFAC,NNB,XI,XID,FIRR,FREG,FD,FDR)
*
*
*       Regular & irregular c.m. force & first derivative.
*       --------------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  XI(3),XID(3),FIRR(3),FREG(3),DV(3),FD(3),FDR(3)
*
*
*       Set non-zero indicator for perturbed c.m.
      NP = 0
      IF (I.GT.N) THEN
          IPAIR = I - N
          IF (LIST(1,2*IPAIR-1).GT.0) NP = 1
      END IF
*
*       Prepare case of single particle or unperturbed c.m. (second call).
      IF (I.LE.N.OR.NP.EQ.0) THEN
*       Copy all KS pairs to JLIST and find MAX(R) for joint treatment.
          NNB1 = NPAIRS
          RMAX1 = 0.0
          DO 10 LJ = 1,NNB1
              JLIST(LJ) = N + LJ
              RMAX1 = MAX(RMAX1,R(LJ))
   10     CONTINUE
*
*       Adopt adequate square distance for c.m. approximation.
          RCM2 = MAX(RCRIT2,CMSEP2*RMAX1**2)
*       Define dummy indices for skipping perturber test.
          JP = 0
          LP = 1
          GO TO 25
      END IF
*
*       Specify variables for treatment of perturbed c.m. particle.
      I2 = 2*IPAIR
      I1 = I2 - 1
      RPERT2 = CMSEP2*R(IPAIR)**2
      BODYIN = 1.0/BODY(I)
*       Initialize perturber list for decision-making.
      NP = LIST(1,I1)
      LP = 2
      JP = LIST(2,I1)
*
*       Use fast force loop for particles satisfying c.m. approximation.
      RCM2 = MAX(RCRIT2,RPERT2)
      NNB1 = 0
      DO 20 J = IFIRST,NTOT
          A1 = X(1,J) - XI(1)
          A2 = X(2,J) - XI(2)
          A3 = X(3,J) - XI(3)
          DV(1) = XDOT(1,J) - XID(1)
          DV(2) = XDOT(2,J) - XID(2)
          DV(3) = XDOT(3,J) - XID(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
*       Form a list of particles for more careful consideration.
          IF (RIJ2.LT.RCM2) THEN
              NNB1 = NNB1 + 1
              JLIST(NNB1) = J
              GO TO 20
          END IF
*
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          FREG(1) = FREG(1) + A1*DR3I
          FREG(2) = FREG(2) + A2*DR3I
          FREG(3) = FREG(3) + A3*DR3I
          FDR(1) = FDR(1) + (DV(1) - A1*DRDV)*DR3I
          FDR(2) = FDR(2) + (DV(2) - A2*DRDV)*DR3I
          FDR(3) = FDR(3) + (DV(3) - A3*DRDV)*DR3I
   20 CONTINUE
*
*       Begin dual purpose force loop (all RIJ2 < RCM2, J > N or I <= N).
   25 DO 60 LJ = 1,NNB1
          JDUM = JLIST(LJ)
          A1 = X(1,JDUM) - XI(1)
          A2 = X(2,JDUM) - XI(2)
          A3 = X(3,JDUM) - XI(3)
          DV(1) = XDOT(1,JDUM) - XID(1)
          DV(2) = XDOT(2,JDUM) - XID(2)
          DV(3) = XDOT(3,JDUM) - XID(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*       First see if the distance exceeds c.m. approximation limit.
          IF (RIJ2.GT.RCM2) GO TO 56
*
          J = JDUM
*       Check whether particle #J satisfies neighbour criteria.
          IF (RIJ2.GT.RCRIT2) GO TO 54
*
*       Consider small step particle (may give large correction terms).
          IF (RIJ2.GT.RS2.AND.STEP(J).GT.SMIN) THEN
              A7 = A1*DV(1) + A2*DV(2) + A3*DV(3)
*       Accept member if maximum penetration factor exceeds 8 per cent.
              IF (A7.GT.VRFAC) GO TO 54
          END IF
          IF (JDUM.EQ.I) GO TO 60
*
*       Obtain force due to current neighbours.
          NNB = NNB + 1
          ILIST(NNB) = J
          KCM = 1
*
*       Advance lower perturber index (includes possible old neighbour).
   26     IF (LP.LE.NP.AND.J.GT.JP) THEN
              LP = LP + 1
              JP = LIST(LP,I1)
*       Include rare case of two consecutive previous neighbours.
              GO TO 26
          END IF
*
*       Decide appropriate expressions from perturber comparison.
          IF (J.NE.JP) THEN
              IF (J.LE.N) GO TO 30
              IF (RIJ2.GT.CMSEP2*R(J-N)**2) GO TO 30
              KDUM = 2*(J - N) - 1
              IF (LIST(1,KDUM).GT.0) THEN
                  K = KDUM
                  J2 = K + 1
                  GO TO 50
              END IF
          ELSE
*       Treat perturbers more carefully.
              IF (LP.LE.NP) THEN
                  LP = LP + 1
                  JP = LIST(LP,I1)
              END IF
              J2 = 0
              IF (J.GT.N) THEN
                  KDUM = 2*(J - N) - 1
                  IF (LIST(1,KDUM).GT.0) THEN
                      J = KDUM
                      J2 = J + 1
                  END IF
              END IF
              GO TO 40
          END IF
*
   30     DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          FIRR(1) = FIRR(1) + A1*DR3I
          FIRR(2) = FIRR(2) + A2*DR3I
          FIRR(3) = FIRR(3) + A3*DR3I
          FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
          FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
          FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
          GO TO 60
*
*       Obtain relevant force on c.m (KCM = 0 denotes regular force).
   40     K = J
   42     L = I1
*       Individual components I1 & I2 are resolved in routine INTGRT.
   45     A1 = X(1,K) - X(1,L)
          A2 = X(2,K) - X(2,L)
          A3 = X(3,K) - X(3,L)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          DV(1) = XDOT(1,K) - XDOT(1,L)
          DV(2) = XDOT(2,K) - XDOT(2,L)
          DV(3) = XDOT(3,K) - XDOT(3,L)
*
          DR2I = 1.0/RIJ2
          DR3I = BODY(K)*BODY(L)*DR2I*SQRT(DR2I)*BODYIN
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          IF (KCM.NE.0) THEN
              FIRR(1) = FIRR(1) + A1*DR3I
              FIRR(2) = FIRR(2) + A2*DR3I
              FIRR(3) = FIRR(3) + A3*DR3I
              FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
              FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
              FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
          ELSE
              FREG(1) = FREG(1) + A1*DR3I
              FREG(2) = FREG(2) + A2*DR3I
              FREG(3) = FREG(3) + A3*DR3I
              FDR(1) = FDR(1) + (DV(1) - A1*DRDV)*DR3I
              FDR(2) = FDR(2) + (DV(2) - A2*DRDV)*DR3I
              FDR(3) = FDR(3) + (DV(3) - A3*DRDV)*DR3I
          END IF
*
          L = L + 1
          IF (L.EQ.I2) GO TO 45
          K = K + 1
          IF (K.EQ.J2) GO TO 42
          GO TO 60
*
*       Treat c.m. approximation for #I and #K as single or composite.
   50     A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          DV(1) = XDOT(1,K) - XID(1)
          DV(2) = XDOT(2,K) - XID(2)
          DV(3) = XDOT(3,K) - XID(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
          DR2I = 1.0/RIJ2
          DR3I = BODY(K)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          IF (KCM.NE.0) THEN
              FIRR(1) = FIRR(1) + A1*DR3I
              FIRR(2) = FIRR(2) + A2*DR3I
              FIRR(3) = FIRR(3) + A3*DR3I
              FD(1) = FD(1) + (DV(1) - A1*DRDV)*DR3I
              FD(2) = FD(2) + (DV(2) - A2*DRDV)*DR3I
              FD(3) = FD(3) + (DV(3) - A3*DRDV)*DR3I
          ELSE
              FREG(1) = FREG(1) + A1*DR3I
              FREG(2) = FREG(2) + A2*DR3I
              FREG(3) = FREG(3) + A3*DR3I
              FDR(1) = FDR(1) + (DV(1) - A1*DRDV)*DR3I
              FDR(2) = FDR(2) + (DV(2) - A2*DRDV)*DR3I
              FDR(3) = FDR(3) + (DV(3) - A3*DRDV)*DR3I
          END IF
*
          K = K + 1
          IF (K.EQ.J2) GO TO 50
          GO TO 60
*
*       Define regular force indicator.
   54 KCM = 0
*       Distinguish between second and first call (I > N & I <= N, J > N)
      IF (JP.EQ.0) THEN
*       Note that first case is for J > N and #I single or unperturbed c.m.
          IF (RIJ2.LT.CMSEP2*R(J-N)**2) THEN
              J2 = 2*(J - N)
              K = J2 - 1
              GO TO 50
          END IF
      ELSE IF (J.LE.N) THEN
*       Consider case of single #J and perturbed c.m.
          IF (RIJ2.LT.RPERT2) THEN
              J2 = 0
              GO TO 40
          END IF
      ELSE
*       Split final case I > N & J > N into two parts according to RPERT2.
          IF (RIJ2.GT.RPERT2) THEN
              IF (RIJ2.GT.CMSEP2*R(J-N)**2) THEN
                  K = J
                  J2 = 0
              ELSE
                  J2 = 2*(J - N)
                  K = J2 - 1
              END IF
*       Adopt c.m. approximation for #I.
              GO TO 50
          ELSE
*       See whether both c.m. bodies should be resolved.
              IF (RIJ2.GT.CMSEP2*R(J-N)**2) THEN
                  J2 = 0
                  GO TO 40
              ELSE
                  J2 = 2*(J - N)
                  K = J2 - 1
                  GO TO 42
              END IF
          END IF
      END IF
*
*       Obtain the regular force due to single body or c.m. particle.
   56     DR2I = 1.0/RIJ2
          DR3I = BODY(JDUM)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          FREG(1) = FREG(1) + A1*DR3I
          FREG(2) = FREG(2) + A2*DR3I
          FREG(3) = FREG(3) + A3*DR3I
          FDR(1) = FDR(1) + (DV(1) - A1*DRDV)*DR3I
          FDR(2) = FDR(2) + (DV(2) - A2*DRDV)*DR3I
          FDR(3) = FDR(3) + (DV(3) - A3*DRDV)*DR3I
   60 CONTINUE
*
*       Check force correction due to regularized chain (same as CMFIRR).
      IF (I.GT.N.AND.NCH.GT.0) THEN
          IF (JP.GT.0) THEN
              CALL KCPERT(I,I1,FIRR,FD)
          END IF
      END IF 
*
      RETURN
*
      END
