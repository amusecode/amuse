      SUBROUTINE CORE
*
*
*       Density centre & core radius.
*       -----------------------------
*
      INCLUDE 'common6.h'
      REAL*8  RLIST(LMAX),RHO,RHOS
      DATA RHOM/1.0D0/
*
*
*       Set generous cutoff.
      RCORE2 = RSCALE**2
      IF (TIME.GT.0.0D0) RCORE2 = MAX(9.0D0*RC**2,RCORE2)
*
*       Select > N/5 central particles.
*       Initialize neighbour densities of all singles and c.m.
    5 NC = 0
      DO 10 I = IFIRST,NTOT
          RHO(I) = 0.D0
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                                   (X(3,I) - RDENS(3))**2
          IF (RI2.LT.RCORE2) THEN
              NC = NC + 1
              JLIST(NC) = I
          END IF
   10 CONTINUE
*

      IF (NC.LT.N/5) THEN
          RCORE2 = 1.5*RCORE2
          GO TO 5
      END IF
*
*       Obtain individual densities.
      RHOS = 0.0
      DO 50 L = 1,NC
          I = JLIST(L)
          NNB = LIST(1,I)
*       Skip and set density to zero of there are no neighbours
          IF (NNB.EQ.0) THEN
          RHO(I) = 0.0D0
          ELSE
*       Estimate D6 to limit the candidates (additional factor of 2).
          D6 = 2.0*(6.0/FLOAT(NNB))**0.66*RS(I)**2
          XI = X(1,I)
          YI = X(2,I)
          ZI = X(3,I)
   20     N6 = 0
*
          DO 25 LJ = 1,NNB
              J = LIST(LJ+1,I)
              RIJ2 = (XI - X(1,J))**2 + (YI - X(2,J))**2 +
     &                                  (ZI - X(3,J))**2
              IF (RIJ2.LT.D6) THEN
                  N6 = N6 + 1
                  ILIST(N6) = J
                  RLIST(N6) = RIJ2
              END IF
   25     CONTINUE
*
          IF (N6.LT.6) THEN
*       Make another iteration in rare cases.
              D6 = 1.5*D6
              IF (N6.LT.NNB) GO TO 20
          END IF
*
*       Sort list of square distances.
          DO 40 II = 1,N6
              DO 35 JJ = II+1,N6
                  IF (RLIST(JJ).LT.RLIST(II)) THEN
                      RDUM = RLIST(II)
                      IDUM = ILIST(II)
                      RLIST(II) = RLIST(JJ)
                      ILIST(II) = ILIST(JJ)
                      RLIST(JJ) = RDUM
                      ILIST(JJ) = IDUM
                  END IF
   35         CONTINUE
   40     CONTINUE
*
          I6 = MIN(N6,6)
*       Calculate total mass of five nearest neighbours.
          XMASS = 0.0D0
          DO 45 LJ = 1,I6-1
              XMASS = XMASS + BODY(ILIST(LJ))
   45     CONTINUE
*
          IF (RLIST(I6).GT.0.0D0) THEN
              RINV32 = 1.0D0/(RLIST(I6)*SQRT(RLIST(I6)))
          ELSE
              RINV32 = 0.0D0
          END IF
          RHO(I) = XMASS*RINV32
*       For multi-mass store also five neighbour particle density
          XNDBL(I) = RINV32
*       Assign zero weight if not enough neighbours.
          IF (N6.LT.6) RHO(I) = 0.0D0
          RHOS = MAX(RHOS,RHO(I))
          END IF
   50 CONTINUE
*
*       Determine the density centre.
      DO 60 K = 1,3
          RDENS(K) = 0.0D0
   60 CONTINUE
*
      RHO1 = 0.0D0
      DO 70 L = 1,NC
          I = JLIST(L)
          DO 65 K = 1,3
              RDENS(K) = RDENS(K) + RHO(I)*X(K,I)
   65     CONTINUE
          RHO1 = RHO1 + RHO(I)
   70 CONTINUE
*
*       Set current density centre based on improved determination.
      RHO1 = MAX(RHO1,ZMASS/RSCALE**3)
      DO 75 K = 1,3
          RDENS(K) = RDENS(K)/RHO1
   75 CONTINUE
*
*       Obtain density radius & averaged density.
      RC = 0.0D0
      RHO2 = 0.0D0
      DO 80 L = 1,NC
          I = JLIST(L)
          XID = X(1,I) - RDENS(1)
          YID = X(2,I) - RDENS(2)
          ZID = X(3,I) - RDENS(3)
          RID2 = XID**2 + YID**2 + ZID**2
          RC = RC + RHO(I)**2*RID2
          RHO2 = RHO2 + RHO(I)**2
   80 CONTINUE
*
*       Form core radius and average & maximum density (scaled in OUTPUT).
      IF (RHO2.GT.0.0D0) RC = SQRT(RC/RHO2)
      RHOD = RHO2/RHO1
      RHOD = (3.0/12.566)*RHOD
      RHOM = (3.0/12.566)*RHOS
*       NB! Scaling factor 5 of Casertano & Hut eq. (V.3) replaced by XMASS.
      IF (RC.EQ.0.0D0) THEN
          RC = RSCALE
          NC = N/2
      END IF
*
*       Set square core radius and its inverse for routine REGINT.
      RC2 = RC**2
      RC2IN = 1.0/RC2
*
*       Sum particles & mass inside the core radius and set rms velocity.
      NC1 = 0
      VC = 0.D00
      ZMC = 0.D00
      DO 85 L = 1,NC
          I = JLIST(L)
          XID = X(1,I) - RDENS(1)
          YID = X(2,I) - RDENS(2)
          ZID = X(3,I) - RDENS(3)
          RID2 = XID**2 + YID**2 + ZID**2
          IF (RID2.LT.RC2) THEN
              NC1 = NC1 + 1
              VC = VC + XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
              ZMC = ZMC + BODY(I)
          END IF
   85 CONTINUE
*
*       Set core membership & rms velocity.
      NC = MAX(NC1,2)
      VC = SQRT(VC/FLOAT(NC))
*
      RETURN
*
      END
