      SUBROUTINE CSTAB3(ITERM)
*
*
*       Perturbed triple chain stability test.
*       --------------------------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      REAL*8  M,MB,MB1,R2(NMX,NMX),XCM(3),VCM(3),XX(3,3),VV(3,3)
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
      I4 = IJ(4)
      MB = M(I1) + M(I2)
      MB1 = MB + M(I3)
      MB4 = MB1 + M(I4)
*
*       Form output diagnostics.
      VREL2 = 0.0D0
      VREL21 = 0.0D0
      RDOT = 0.0D0
      RDOT3 = 0.0D0
      RI2 = 0.0D0
      R4 = 0.0D0
      V4 = 0.0
      RDOT4 = 0.0D0
      DO 10 K = 1,3
          J1 = 3*(I1-1) + K
          J2 = 3*(I2-1) + K
          J3 = 3*(I3-1) + K
          J4 = 3*(I4-1) + K
          VREL2 = VREL2 + (V(J1) - V(J2))**2
          RDOT = RDOT + (X(J1) - X(J2))*(V(J1) - V(J2))
          XCM(K) = (M(I1)*X(J1) + M(I2)*X(J2))/MB
          VCM(K) = (M(I1)*V(J1) + M(I2)*V(J2))/MB
          RI2 = RI2 + (X(J3) - XCM(K))**2
          VREL21 = VREL21 + (V(J3) - VCM(K))**2
          RDOT3 = RDOT3 + (X(J3) - XCM(K))*(V(J3) - VCM(K))
          XC3 = (M(I3)*X(J3) + MB*XCM(K))/MB1
          VC3 = (M(I3)*V(J3) + MB*VCM(K))/MB1
          R4 = R4 + (X(J4) - XC3)**2
          V4 = V4 + (V(J4) - VC3)**2
          RDOT4 = RDOT4 + (X(J4) - XC3)*(V(J4) - VC3)
          XX(K,1) = X(J1)
          XX(K,2) = X(J2)
          XX(K,3) = X(J3)
          VV(K,1) = V(J1)
          VV(K,2) = V(J2)
          VV(K,3) = V(J3)
   10 CONTINUE
*
*       Evaluate orbital elements for inner and outer motion.
      RB = SQRT(R2(I1,I2))
      R3 = SQRT(RI2)
      SEMI = 2.0D0/RB - VREL2/MB
      SEMI = 1.0/SEMI
      ECC = SQRT((1.0D0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB))
      SEMI1 = 2.0/R3 - VREL21/MB1
      SEMI1 = 1.0/SEMI1
      ECC1 = SQRT((1.0D0 - R3/SEMI1)**2 + RDOT3**2/(SEMI1*MB1))
      R4 = SQRT(R4)
      SEMI4 = 2.0/R4 - V4/MB4
      SEMI4 = 1.0/SEMI4
      ECC4 = SQRT((1.0 - R4/SEMI4)**2 + RDOT4**2/(SEMI4*MB4))
*
*       Form hierarchical stability ratio (Eggleton & Kiseleva 1995).
      QL = MB/M(I3)
      Q1 = MAX(M(I2)/M(I1),M(I1)/M(I2))
      Q3 = QL**0.33333
      Q13 = Q1**0.33333
      AR = 1.0 + 3.7/Q3 - 2.2/(1.0 + Q3) + 1.4/Q13*(Q3 - 1.0)/(Q3 + 1.0)
*
      EK = AR*SEMI*(1.0D0 + ECC)
      PMIN = SEMI1*(1.0D0 - ECC1)
*
*       Replace the EK criterion by the MA analytical stability formula.
      Q0 = M(I3)/MB
      IF (ECC1.LT.1.0) THEN
          XFAC = (1.0 + Q0)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
      ELSE
          XFAC = 40.0*(1.0 + Q0)
      END IF
      FE = 1.0
      PCRIT = 2.8*FE*XFAC**0.4*SEMI
*
*       Obtain the inclination.
      CALL INCLIN(XX,VV,XCM,VCM,ALPHA)
*
*       Include fudge factor for inclination effect.
      YFAC = 1.0 - 0.3*ALPHA/3.14
      PCRIT = YFAC*PCRIT
*
*       Check hierarchical stability condition for 3 + 1 configuration.
      ITERM = 0
      IF (PMIN.GT.PCRIT.AND.SEMI1.GT.0.0.AND.SEMI4.GT.0.0.AND.
     &    RB.GT.SEMI) THEN
          G4 = 2.0*M(I4)/MB4*(R3/R4)**3
          PMIN4 = SEMI4*(1.0 - ECC4)
          ZFAC = (1.0 + M(I4)/MB1)*(1.0 + ECC4)/SQRT(1.0 - ECC4)
          ZCRIT = 2.8*FE*ZFAC**0.4*SEMI4
          PZ = PMIN4/ZCRIT
*
*       Ignore outermost body for perturbation test of inner triple.
          IF (N.EQ.5) THEN
              G4 = 2.0*M(I4)/MB1*(R3/PMIN4)**3
              IF (G4.LT.0.1) THEN
                  if(rank.eq.0)
     &            WRITE (6,19)  NSTEP1, ECC4, SEMI, PCRIT, PMIN4, G4
   19             FORMAT (' CSTAB3    # E4 A PCR PM4 G4 ',
     &                                I6,F8.4,1P,4E10.2)
                  ITERM = -1
              END IF
              GO TO 40
          END IF
*
*       Continue chain integration if outer orbit unstable or large pert.
          IF (PMIN4.LT.ZCRIT.OR.G4.GT.0.2) GO TO 40
*
          ITERM = -1
          if(rank.eq.0)
     &    WRITE (6,20)  ECC, ECC1, SEMI, SEMI1, PMIN, PCRIT, EK, PZ
   20     FORMAT (' CSTAB3    E =',F6.3,'  E1 =',F6.3,'  A =',1P,E8.1,
     &                     '  A1 =',E8.1,'  PM =',E9.2,'  PC =',E9.2,
     &                     '  EK =',E9.2,'  P/Z =',0P,F5.2)
          RI = SQRT(CM(1)**2 + CM(2)**2 + CM(3)**2)
          EMAX = 0.0
          if(rank.eq.0)
     &    WRITE (81,30)  TIMEC, RI, NAMEC(I3), QL, Q1, ECC, ECC1,
     &                   SEMI, SEMI1, PCRIT/PMIN, 180.*ALPHA/3.14, EMAX
   30     FORMAT (F9.5,F5.1,I6,2F6.2,2F6.3,1P,2E10.2,0P,F5.2,F6.1,F6.3)
          CALL FLUSH(81)
      END IF
*
   40 RETURN
*
      END
