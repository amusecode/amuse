      SUBROUTINE CHSTAB(ITERM)
*
*
*       Chain stability test.
*       ---------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      REAL*8  M,MB,MB1,R2(NMX,NMX),XCM(3),VCM(3),XX(3,3),VV(3,3),
     &        A1(3),A2(3),XREL(3),VREL(3),EI(3),HI(3),HO(3)
      INTEGER  IJ(NMX)
      INCLUDE "mpif.h"
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
*
*
*       Sort particle separations (I1 & I2 form closest pair).
      CALL R2SORT(IJ,R2)
      I1 = IJ(1)
      I2 = IJ(2)
      I3 = IJ(3)
      MB = M(I1) + M(I2)
      MB1 = MB + M(I3)
*
*       Form output diagnostics.
      VREL2 = 0.0D0
      VREL21 = 0.0D0
      RDOT = 0.0D0
      RDOT3 = 0.0D0
      R3 = 0.0D0
      DO 5 K = 1,3
          J1 = 3*(I1-1) + K
          J2 = 3*(I2-1) + K
          J3 = 3*(I3-1) + K
          VREL2 = VREL2 + (V(J1) - V(J2))**2
          RDOT = RDOT + (X(J1) - X(J2))*(V(J1) - V(J2))
          XCM(K) = (M(I1)*X(J1) + M(I2)*X(J2))/MB
          VCM(K) = (M(I1)*V(J1) + M(I2)*V(J2))/MB
          R3 = R3 + (X(J3) - XCM(K))**2
          VREL21 = VREL21 + (V(J3) - VCM(K))**2
          RDOT3 = RDOT3 + (X(J3) - XCM(K))*(V(J3) - VCM(K))
          XX(K,1) = X(J1)
          XX(K,2) = X(J2)
          XX(K,3) = X(J3)
          VV(K,1) = V(J1)
          VV(K,2) = V(J2)
          VV(K,3) = V(J3)
    5 CONTINUE
*
*       Evaluate orbital elements for inner and outer motion.
      RB = SQRT(R2(I1,I2))
      R3 = SQRT(R3)
      SEMI = 2.0D0/RB - VREL2/MB
      SEMI = 1.0/SEMI
      ECC = SQRT((1.0D0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB))
      SEMI1 = 2.0/R3 - VREL21/MB1
      SEMI1 = 1.0/SEMI1
      ECC1 = SQRT((1.0D0 - R3/SEMI1)**2 + RDOT3**2/(SEMI1*MB1))
      PMIN = SEMI1*(1.0D0 - ECC1)
*
*       Form hierarchical stability ratio (Eggleton & Kiseleva 1995).
*     QL = MB/M(I3)
*     Q1 = MAX(M(I2)/M(I1),M(I1)/M(I2))
*     Q3 = QL**0.33333
*     Q13 = Q1**0.33333
*     AR = 1.0 + 3.7/Q3 - 2.2/(1.0 + Q3) + 1.4/Q13*(Q3 - 1.0)/(Q3 + 1.0)
*     EK = AR*SEMI*(1.0D0 + ECC)
*
*       Obtain the inclination.
      CALL INCLIN(XX,VV,XCM,VCM,ALPHA)
*
*       Replace the EK criterion by the MA 1999 stability expression.
*     PCRIT = stability(M(I1),M(I2),M(I3),ECC,ECC1,ALPHA)*SEMI
*       Add 1% for perturbation to avoid repeated switching.
*     PCRIT = 1.01*PCRIT
*
*       Evaluate the general stability function (Mardling, Cambody 2008).
      IF (ECC1.LT.1.0.AND.ECC.LT.1.0) THEN
          NST = NSTAB(SEMI,SEMI1,ECC,ECC1,ALPHA,M(I1),M(I2),M(I3))
          IF (NST.EQ.0) THEN
              PCRIT = 0.99*PMIN
          ELSE
              PCRIT = 1.01*PMIN
          END IF
      ELSE
*       Set nominal failed value for SEMI1 > 0 test.
          PCRIT = 1.01*PMIN
      END IF
*
*       Prepare evaluation of maximum eccentricity (see routine HIGROW).
      DO 10 K = 1,3
          XREL(K) = XX(K,1) - XX(K,2)
          VREL(K) = VV(K,1) - VV(K,2)
   10 CONTINUE
      A12 = 0.0
      A22 = 0.0
      A1A2 = 0.0
      RI2 = 0.0
      VI2 = 0.0
      RVI = 0.0
      DO 12 K = 1,3
          K1 = K + 1
          IF (K1.GT.3) K1 = 1
          K2 = K1 + 1
          IF (K2.GT.3) K2 = 1
          A1(K) = XREL(K1)*VREL(K2) - XREL(K2)*VREL(K1)
          A2(K) = (XX(K1,3) - XCM(K1))*(VV(K2,3) - VCM(K2))
     &          - (XX(K2,3) - XCM(K2))*(VV(K1,3) - VCM(K1))
          A12 = A12 + A1(K)**2
          A22 = A22 + A2(K)**2
          A1A2 = A1A2 + A1(K)*A2(K)
          RI2 = RI2 + XREL(K)**2
          VI2 = VI2 + VREL(K)**2
          RVI = RVI + XREL(K)*VREL(K)
   12 CONTINUE
*
*       Construct the Runge-Lenz vector (Heggie & Rasio 1995, Eq.(5)).
      EI2 = 0.0
      DO 15 K = 1,3
          EI(K) = (VI2*XREL(K) - RVI*VREL(K))/MB - XREL(K)/SQRT(RI2)
          EI2 = EI2 + EI(K)**2
   15 CONTINUE
      EI2 = MIN(EI2,0.9999D0)
*
*       Define unit vectors for inner eccentricity and angular momenta.
      COSJ = 0.0
      SJSG = 0.0
      DO 18 K = 1,3
          EI(K) = EI(K)/SQRT(EI2)
          HI(K) = A1(K)/SQRT(A12)
          HO(K) = A2(K)/SQRT(A22)
          COSJ = COSJ + HI(K)*HO(K)
          SJSG = SJSG + EI(K)*HO(K)
   18 CONTINUE
*
*       Evaluate the expressions A & Z.
      A = COSJ*SQRT(1.0 - EI2)
      Z = (1.0 - EI2)*(2.0 - COSJ**2) + 5.0*EI2*SJSG**2
*
*       Obtain maximum inner eccentricity (Douglas Heggie, Sept. 1995).
      Z2 = Z**2 + 25.0 + 16.0*A**4 - 10.0*Z - 20.0*A**2 - 8.0*A**2*Z
      Z2 = MAX(Z2,0.0D0)
      EMAX = (Z + 1.0 - 4.0*A**2 + SQRT(Z2))/6.0
      EMAX = MAX(EMAX,0.0001D0)
      EMAX = SQRT(EMAX)
*
*       Check hierarchical stability condition for bound close pair (RB > a).
      ITERM = 0
      ALPHA = 180.0*ALPHA/3.1415
      IF (PMIN.GT.PCRIT.AND.SEMI.GT.0.0.AND.SEMI1.GT.0.0.AND.
     &    RB.GT.SEMI) THEN
          IF (RDOT3.GT.0.0.AND.R3.GT.3.0*SEMI*(1.0 + ECC1)) THEN
          ITERM = -1
          if(rank.eq.0)
     &    WRITE (6,20)  NAMEC(I1), NAMEC(I2), NAMEC(I3), ECC, EMAX,
     &                  ECC1, SEMI, SEMI1, PMIN, PCRIT, ALPHA
   20     FORMAT (' NEW HIARCH    NM =',3I6,'  E =',F6.3,'  EX =',F7.4,
     &                         '  E1 =',F6.3,'  A =',1P,E8.1,
     &                         '  A1 =',E8.1,'  PM =',E9.2,
     &                         '  PC =',E9.2,'  IN =',0P,F7.1)
          RI = SQRT(CM(1)**2 + CM(2)**2 + CM(3)**2)
          Q0 = M(I3)/MB
          if(rank.eq.0)
     &    WRITE (81,30)  TIMEC, RI, NAMEC(I3), Q0, ECC, EMAX, ECC1,
     &                   SEMI, SEMI1, PCRIT/PMIN, ALPHA
   30     FORMAT (F8.1,F5.1,I6,F6.2,3F6.3,1P,2E10.2,0P,F5.2,F8.1)
          CALL FLUSH(81)
          END IF
*       Include termination test for wide triple system (exclude ECC1 > 0.9).
      ELSE IF (PMIN.GT.3.0*SEMI*(1.0 + ECC).AND.SEMI.GT.0.0.AND.
     &         ECC1.LT.0.9) THEN
*       Wait for favourable configuration (R > SEMI, R3 > 3*SEMI & RDOT3 > 0.
          IF (RB.GT.SEMI.AND.R3.GT.SEMI1.AND.RDOT3.GT.0.0) THEN
              APO = SEMI*(1.0 + ECC)
              if(rank.eq.0)
     &        WRITE (6,40)  ECC, ECC1, ALPHA, RB, R3, PCRIT, PMIN, APO
   40         FORMAT (' WIDE CHAIN    E E1 IN RB R3 PC PM APO ',
     &                                2F7.3,F7.1,1P,5E10.2)
              CALL FLUSH(6)
              ITERM = -1
          END IF
*       Enforce termination of long-lived configuration using two limits.
      ELSE IF (NSTEP1.GT.50000.AND.PMIN.GT.0.9*PCRIT.AND.
     &    R3.GT.SEMI1) THEN
          ITERM = -1
      ELSE IF (NSTEP1.GT.500000.AND.R3.GT.SEMI1) THEN
          ITERM = -1
      END IF
*
      RETURN
*
      END
