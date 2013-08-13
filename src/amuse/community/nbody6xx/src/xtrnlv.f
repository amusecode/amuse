      SUBROUTINE XTRNLV(I1,I2)
*
*
*       External potential and virial energy.
*       -------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2
      REAL*8 XI(3),XIDOT(3),FIRR(3),FREG(3),FD(3),FDR(3)
      SAVE FIRST
      LOGICAL FIRST
      DATA FIRST /.TRUE./
*
*
      ET = 0.0D0
*       Skip external potential during scaling (parameters not defined).
      IF (FIRST.AND.TIME.EQ.0.0D0) THEN
          FIRST = .FALSE.
          GO TO 30
*       Treat all cases uniformly (but note Coliolis term in FIRR).
      ELSE IF (KZ(14).LE.4.AND.I2.NE.I1) THEN
          VIR = 0.0
          DO 5 I = I1,I2
              DO 2 K = 1,3
*       Initialize scalars for force contributions.
                  XI(K) = X(K,I)
                  XIDOT(K) = XDOT(K,I)
                  FIRR(K) = 0.0
                  FREG(K) = 0.0
                  FD(K) = 0.0
                  FDR(K) = 0.0
    2         CONTINUE
*       Evaluate the virial energy using definition sum {m*r*F}.
              CALL XTRNLF(XI,XIDOT,FIRR,FREG,FD,FDR,1)
*       Note that FIRR(K) only contributes for #14 = 1 or 2.
              DO 4 K = 1,3
                  VIR = VIR + BODY(I)*XI(K)*(FIRR(K) + FREG(K))
    4         CONTINUE
    5     CONTINUE
*       Quit for pure 3D tidal field (case of full N summation).
          IF (KZ(14).EQ.3) GO TO 40
      END IF
*
*       See whether to include a linearized galactic tidal force.
      IF (KZ(14).LE.2) THEN
          DO 10 I = I1,I2
              ET = ET - 0.5D0*BODY(I)*(TIDAL(1)*X(1,I)**2 +
     &                                 TIDAL(3)*X(3,I)**2)
   10     CONTINUE
      ELSE IF (KZ(14).EQ.3) THEN
*       Retain 3D potential energy expressions for future use.
          I = I1
          RG2 = 0.0
          RI2 = 0.0
          DO 20 K = 1,3
              XI(K) = X(K,I)
              RG2 = RG2 + RG(K)**2
              RI2 = RI2 + (RG(K) + XI(K))**2
   20     CONTINUE
*
*       Include galaxy point mass term for body #I in differential form.
          IF (GMG.GT.0.0D0) THEN
              ET = ET + GMG*(1.0/SQRT(RI2) - 1.0/SQRT(RG2))
          END IF
*
*       Add optional Miyamoto disk potential.
          IF (DISK.GT.0.0D0) THEN
              R2 = (RG(1) + XI(1))**2 + (RG(2) + XI(2))**2
              BZ = SQRT(B**2 + (RG(3) + XI(3))**2)
              AZ = SQRT(R2 + (A + BZ)**2)
              R20 = RG(1)**2 + RG(2)**2
              BZ0 = SQRT(B**2 + RG(3)**2)
              AZ0 = SQRT(R20 + (A + BZ0)**2)
              ET = ET + DISK*(1.0/AZ - 1.0/AZ0)
          END IF
*
*       Check addition of differential logarithmic potential.
          IF (V02.GT.0.0D0) THEN
              ET = ET + 0.5*V02*(LOG(RI2) - LOG(RG2))
          END IF
*       Form the differential potential energy due to tides.
          ET = BODY(I)*ET
      END IF
*
*       Place sum in ETIDE and single particle contribution in HT.
   30 IF (I2.GT.I1) THEN
          ETIDE = ET
      ELSE
          HT = ET
      END IF
*
   40 RETURN
*
      END
