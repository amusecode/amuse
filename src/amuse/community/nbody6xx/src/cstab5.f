      SUBROUTINE CSTAB5(ITERM)
*
*
*       Five-body chain stability test.
*       -------------------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      REAL*8  M,MB,MB1,MB2,R2(NMX,NMX),XCM(3),VCM(3),XX(3,3),VV(3,3),
     &        XCM2(3),VCM2(3),XB(3),VB(3),M1,M2,M3
      INTEGER  IJ(NMX),ISORT(NMX)
      INCLUDE "mpif.h"
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
*
*
*       Exclude unlikely stability configuration (ratio 1/4 - 4).
      IR = 0
      DO 1 I = 1,N-2
          RATIO = RINV(I+1)/RINV(I)
          IF (RATIO.GT.0.25.AND.RATIO.LT.4.0) THEN
              IR = IR + 1
          END IF
    1 CONTINUE
*
*       Exit on at least one unfavourable distance ratio.
      IF (IR.GT.0) THEN
*         if(rank.eq.0)
*    &    WRITE (6,2)  I, (1.0/RINV(K),K=1,N-1)
*   2     FORMAT (' CSTAB5 TEST    I R ',I4,1P,5E9.1)
          GO TO 40
      END IF
*
*       Ensure a close binary at beginning or end of chain.
*     RB = 0.04*RSUM
*     IF (1.0/RINV(1).GT.RB.AND.1.0/RINV(N-1).GT.RB) GO TO 40
*
*       Check stability of three innermost bodies and one perturber.
      CALL CSTAB3(ITERM)
      IF (ITERM.LT.0) GO TO 40
*
      CALL HPSORT(N-1,RINV,ISORT)
      IB2 = ISORT(2)
      IB3 = ISORT(3)
      RATIO = RINV(IB3)/RINV(IB2)
      IF (RATIO.GT.0.25.AND.RATIO.LT.4.0) GO TO 40
*
*       Sort particle separations (I1 & I2 form closest pair).
      CALL R2SORT(IJ,R2)
*
*       Save indices with I1 & I2 and I3 & I4 as inner & outer binary.
      I1 = IJ(1)
      I2 = IJ(2)
      I3 = IJ(3)
      I4 = IJ(4)
      I5 = 0
      IF (N.EQ.2) THEN
          I3 = I2
          I4 = I1
      ELSE IF (N.EQ.3) THEN
          I4 = I1
*       Note N = 4 is already defined correctly (but redundant here).
      ELSE IF (N.GT.4) THEN
*       Determine indices of second closest pair (avoid pair I1-I2).
          RX1 = 1.0
          RX0 = R2(I1,I2)
          DO 5 J1 = 1,N
              IF (J1.EQ.I1.OR.J1.EQ.I2) GO TO 5
              DO 4 J2 = J1+1,N
                  IF (J2.EQ.I1.OR.J2.EQ.I2) GO TO 4
                  IF (R2(J1,J2).LT.RX1.AND.R2(J1,J2).GT.RX0) THEN
                      RX1 = R2(J1,J2)
                      I3 = J1
                      I4 = J2
                  END IF
    4         CONTINUE
    5     CONTINUE
*       Identify remaining single particle(s) by exclusion.
          DO 8 I = 1,N
              IF (I.EQ.I1.OR.I.EQ.I2.OR.I.EQ.I3.OR.I.EQ.I4) GO TO 8
              IF (I5.EQ.0) THEN
                  I5 = I
              ELSE 
                  I6 = I
              END IF
    8     CONTINUE
      END IF
*
      MB = M(I1) + M(I2)
      MB2 = M(I3) + M(I4)
      M5 = M(I5)
      MB1 = MB + MB2 + M5
*
*       Form orbital parameters with c.m. of I3 & I4 as third body.
      VREL2 = 0.0D0
      VREL3 = 0.0D0
      RDOT = 0.0D0
      D12 = 0.0
      D15 = 0.0
      D25 = 0.0
      V12 = 0.0
      V15 = 0.0
      V25 = 0.0
      DO 10 K = 1,3
          J1 = 3*(I1-1) + K
          J2 = 3*(I2-1) + K
          J3 = 3*(I3-1) + K
          J4 = 3*(I4-1) + K
          J5 = 3*(I5-1) + K
          VREL2 = VREL2 + (V(J1) - V(J2))**2
          RDOT = RDOT + (X(J1) - X(J2))*(V(J1) - V(J2))
          XCM(K) = (M(I1)*X(J1) + M(I2)*X(J2))/MB
          VCM(K) = (M(I1)*V(J1) + M(I2)*V(J2))/MB
          XCM2(K) = (M(I3)*X(J3) + M(I4)*X(J4))/MB2
          VCM2(K) = (M(I3)*V(J3) + M(I4)*V(J4))/MB2
          VREL3 = VREL3 + (V(J3) - V(J4))**2
          D12 = D12 + (XCM(K) - XCM2(K))**2
          D15 = D15 + (XCM(K) - X(J5))**2
          D25 = D25 + (XCM2(K) - X(J5))**2
          V12 = V12 + (VCM(K) - VCM2(K))**2
          V15 = V15 + (VCM(K) - V(J5))**2
          V25 = V25 + (VCM2(K) - V(J5))**2
   10 CONTINUE
*
*       Select inner binary and form relevant quantities.
      RB1 = 0.0
      VB1 = 0.0
      RD1 = 0.0
      IF (D12.LT.MIN(D15,D25)) THEN
          SEMI = 2.0/SQRT(D12) - V12/(MB + MB2)
*         Q0 = M5/(MB + MB2)
          M1 = MB
          M2 = MB2
          M3 = M5
          DO 11 K = 1,3
              J5 = 3*(I5-1) + K
              XB(K) = (MB*XCM(K) + MB2*XCM2(K))/(MB + MB2)
              VB(K) = (MB*VCM(K) + MB2*VCM2(K))/(MB + MB2)
              RD1 = RD1 + (XB(K) - X(J5))*(VB(K) - V(J5))
              RB1 = RB1 + (XB(K) - X(J5))**2
              VB1 = VB1 + (VB(K) - V(J5))**2
              XX(K,1) = XCM(K)
              VV(K,1) = VCM(K)
              XX(K,2) = XCM2(K)
              VV(K,2) = VCM2(K)
              XX(K,3) = X(J5)
              VV(K,3) = V(J5)
   11     CONTINUE
      ELSE IF (D15.LT.MIN(D12,D25)) THEN
          SEMI = 2.0/SQRT(D15) - V15/(MB + M5)
*         Q0 = MB2/(MB + M5)
          M1 = MB
          M2 = M5
          M3 = MB2
          DO 12 K = 1,3
              J5 = 3*(I5-1) + K
              XB(K) = (MB*XCM(K) + M5*X(J5))/(MB + M5)
              VB(K) = (MB*VCM(K) + M5*V(J5))/(MB + M5)
              RD1 = RD1 + (XB(K) - XCM2(K))*(VB(K) - VCM2(K))
              RB1 = RB1 + (XB(K) - XCM2(K))**2
              VB1 = VB1 + (VB(K) - VCM2(K))**2
              XX(K,1) = XCM(K)
              VV(K,1) = VCM(K)
              XX(K,2) = X(J5)
              VV(K,2) = V(J5)
              XX(K,3) = XCM2(K)
              VV(K,3) = VCM2(K)
   12     CONTINUE
      ELSE
          SEMI = 2.0/SQRT(D25) - V25/(MB2 + M5)
*         Q0 = MB/(MB2 + M5)
          M1 = MB2
          M2 = M5
          M3 = MB
          DO 13 K = 1,3
              J5 = 3*(I5-1) + K
              XB(K) = (MB2*XCM2(K) + M5*X(J5))/(MB2 + M5)
              VB(K) = (MB2*VCM2(K) + M5*V(J5))/(MB2 + M5)
              RD1 = RD1 + (XB(K) - XCM(K))*(VB(K) - VCM(K))
              RB1 = RB1 + (XB(K) - XCM(K))**2
              VB1 = VB1 + (VB(K) - VCM(K))**2
              XX(K,1) = XCM2(K)
              VV(K,1) = VCM2(K)
              XX(K,2) = X(J5)
              VV(K,2) = V(J5)
              XX(K,3) = XCM(K)
              VV(K,3) = VCM(K)
   13     CONTINUE
      END IF
*
*       Specify inner & outer semi-major axis.
      SEMI = 1.0/SEMI
      RB1 = SQRT(RB1)
      SEMI1 = 2.0/RB1 - VB1/MB1
      SEMI1 = 1.0/SEMI1
*
*       Obtain orbital elements for inner and outer motion.
      RB0 = SQRT(R2(I1,I2))
      SEMI0 = 2.0/RB0 - VREL2/MB
      SEMI0 = 1.0/SEMI0
      ECC0 = SQRT((1.0D0 - RB0/SEMI0)**2 + RDOT**2/(SEMI0*MB))
      ECC1 = SQRT((1.0D0 - RB1/SEMI1)**2 + RD1**2/(SEMI1*MB1))
      RB2 = SQRT(R2(I3,I4))
      SEMI2 = 2.0/RB2 - VREL3/MB2
      SEMI2 = 1.0/SEMI2
*
*       Obtain the inclination.
      CALL INCLIN(XX,VV,XB,VB,ALPHA)
*
*       Employ the general stability criterion.
      PCRIT = stability(M1,M2,M3,ECC0,ECC1,ALPHA)*SEMI
*
*       Modify correction factor by widest binary (cf. routine IMPACT).
      IF (SEMI2.GT.0.0) PCRIT = (1.0 + 0.1*SEMI2/SEMI)*PCRIT
*
*       Check hierarchical stability condition for bound close pair.
      PMIN = SEMI1*(1.0D0 - ECC1)
      ITERM = 0
      IF (PMIN.GT.PCRIT.AND.SEMI.GT.0.0.AND.SEMI1.GT.0.0.AND.
     &    RB0.GT.SEMI0.AND.SEMI2.GT.0.0) THEN
          ITERM = -1
          if(rank.eq.0)
     &    WRITE (6,20)  ECC0, ECC1, SEMI, SEMI1, PMIN, PCRIT, SEMI2,
     &                  180.0*ALPHA/3.14
   20     FORMAT (' CSTAB5    E0 =',F6.3,'  E1 =',F6.3,'  A =',1P,E8.1,
     &                     '  A1 =',E8.1,'  PM =',E9.2,'  PC =',E9.2,
     &                     '  A2 =',E8.1,'  IN =',0P,F6.1)
      END IF
*
   40 RETURN
*
      END
