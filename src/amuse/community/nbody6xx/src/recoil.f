      SUBROUTINE RECOIL(IEND)
*
*
*       Binary analysis and final output.
*       ---------------------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      REAL*8  R2(NMX,NMX),RC(3),VC(3),XX(3,3),VV(3,3)
      INTEGER  IJ(NMX)
      INCLUDE "mpif.h"
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/CALLS/  TPR,TKPR,STEP,IDER,ICALL,NFN,NREG,ITER,IMCIRC
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/CLUMP/   BODYS(NMX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NMX,5),ISYS(5)
      SAVE
*
*
*       Sort particle separations (I1 & I2 form closest pair).
      CALL R2SORT(IJ,R2)
      I1 = IJ(1)
      I2 = IJ(2)
      I3 = IJ(3)
      I4 = IJ(4)
      IF (N.EQ.2) THEN
          I3 = I2
          I4 = I1
      ELSE IF (N.EQ.3) THEN
          I4 = I1
      ELSE IF (N.GE.4) THEN
*       Determine indices of second closest pair (avoid pair I1-I2).
          RX1 = 1.0
          RX0 = R2(I1,I2)
          DO 2 J1 = 1,N
              IF (J1.EQ.I1.OR.J1.EQ.I2) GO TO 2
              DO 1 J2 = J1+1,N
                  IF (J2.EQ.I1.OR.J2.EQ.I2) GO TO 1
                  IF (R2(J1,J2).LT.RX1.AND.R2(J1,J2).GT.RX0) THEN
                      RX1 = R2(J1,J2)
                      I3 = J1
                      I4 = J2
                  END IF
    1         CONTINUE
    2     CONTINUE
      END IF
      K1 = I3
      K2 = I4
*
*       Ensure that original KS pair is considered (H > 0 is possible).
      IF (N.LE.3.AND.IEND.EQ.0) THEN
          DO 5 K = 1,N
              IF (INAME(K).EQ.1) K1 = K
              IF (INAME(K).EQ.2) K2 = K
    5     CONTINUE
      END IF
*
*       Form output diagnostics.
      VREL2 = 0.0D0
      VREL21 = 0.0D0
      VREL23 = 0.0D0
      RDOT = 0.0D0
      RC2 = 0.0D0
      VC2 = 0.0D0
      RCDOT = 0.0D0
      MB = M(I1) + M(I2)
*
      DO 10 K = 1,3
          J1 = 3*(I1-1) + K
          J2 = 3*(I2-1) + K
          J3 = 3*(I3-1) + K
          J4 = 3*(I4-1) + K
          L1 = 3*(K1-1) + K
          L2 = 3*(K2-1) + K
          VREL2 = VREL2 + (V(J1) - V(J2))**2
          VREL21 = VREL21 + (V(J3) - V(J4))**2
          VREL23 = VREL23 + (V(L1) - V(L2))**2
          RDOT = RDOT + (X(J1) - X(J2))*(V(J1) - V(J2))
          RC(K) = (M(I1)*X(J1) + M(I2)*X(J2))/MB
          VC(K) = (M(I1)*V(J1) + M(I2)*V(J2))/MB
          RC2 = RC2 + (RC(K) - X(J3))**2
          VC2 = VC2 + (VC(K) - V(J3))**2
          RCDOT = RCDOT + (RC(K) - X(J3))*(VC(K) - V(J3))
          XX(K,1) = X(J1)
          XX(K,2) = X(J2)
          XX(K,3) = X(J3)
          VV(K,1) = V(J1)
          VV(K,2) = V(J2)
          VV(K,3) = V(J3)
   10 CONTINUE
*
*       Evaluate orbital elements.
      RB = SQRT(R2(I1,I2))
      RB1 = SQRT(R2(I3,I4))
      R13 = SQRT(R2(I1,I3))
      R24 = SQRT(R2(I2,I4))
      SEMI = 2.0D0/RB - VREL2/MB
      SEMI = 1.0/SEMI
      ECC = SQRT((1.0D0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB))
      EB = -M(I1)*M(I2)/(2.0D0*SEMI)
      ET = -M(I1)*M(I2)/(2.0D0*SEMI*CM(8))
*
*       Determine parameters for possible hierarchy of body #I3.
      RCP = SQRT(RC2)
      A1 = 2.0/RCP - VC2/(MB + M(I3))
      A1 = 1.0/A1
      ECC1 = SQRT((1.0 - RCP/A1)**2 + RCDOT**2/(A1*(MB + M(I3))))
      PMIN = A1*(1.0 - ECC1)
*
*       Obtain binding energy of K1 & K2 (in case I1 & I2 not bound).
      SEMI2 = 2.0D0/SQRT(R2(K1,K2)) - VREL23/(M(K1) + M(K2))
      SEMI2 = 1.0/SEMI2
      EB2 = -M(K1)*M(K2)/(2.0D0*SEMI2)
*
*       Include possible second binary (suppress positive energy).
      IF (N.GT.3) THEN
          MB1 = M(I3) + M(I4)
          SEMI1 = 2.0D0/RB1 - VREL21/MB1
          SEMI1 = 1.0/SEMI1
          EB1 = -M(I3)*M(I4)/(2.0D0*SEMI1)
          EB1 = MIN(EB1,0.0D0)
      ELSE
          R24 = 1.0E+10
          IF (TIMEC.LE.0.0D0) EB1 = 0.0
      END IF
*
*       Ensure that perturbed boundary exceeds system size.
      RMAX = MAX(SQRT(R2(I1,I3)),SQRT(R2(I2,I4)))
      RMAXC = MAX(RMAXC,1.2*RMAX)
*
*       Save initial energies (note: K1 & K2 may be most energetic).
      IF (IEND.EQ.0) THEN
          E0 = ENERGY
          DE = 0.0
*         IF (TIMEC.GT.0.0D0) GO TO 50
          EBCH0 = CM(9)
          EB0 = MIN(EB,EB2)
          EB10 = EB1
          A0 = SEMI
          ECC0 = ECC
          IF (EB2.LT.EB) A0 = SEMI2
*       Note closest pair (#I1/I2 may be intruder or wide binary at peri).
          IF (N.EQ.3) THEN
              NAME1 = NAMEC(1)
              NAME2 = NAMEC(2)
          ELSE
              NAME1 = NAMEC(I1)
              NAME2 = NAMEC(I2)
          END IF
          GO TO 50
      ELSE IF (IEND.EQ.1) THEN
          DE = DE + ABS(ENERGY - E0)
          E0 = ENERGY
*       Initialize second binary on subsequent absorption.
          IF (EB10.EQ.0.0D0) EB10 = EB1
          GO TO 50
      END IF
*
*       Determine nearest and most distant binary perturber (N = 4).
      IF (R13.LT.R24) THEN
          IM = I3
          RM = R13
          IMAX = I4
          RMAX = R24
      ELSE
          IM = I4
          RM = R24
          IMAX = I3
          RMAX = R13
      END IF
*
*       Estimate perturbations on smallest binary.
      IF (N.GT.3) THEN
          IF (N.GT.4) RM = MAX(R13,R24)
          GB = 2.0*MB1/MB*(RB/RM)**3
          IF (EB1.LT.0.0) GB = GB*M(IM)/MB1
          G4 = 2.0*M(IMAX)/(MB + M(IM))*(RM/RMAX)**3
      ELSE
          GB = M(I3)/MB*(RB/R13)**3
          G4 = 0.0
          EB1 = 0.0
      END IF
*
*       Form net binary energy change (added to CHCOLL in routine CHTERM).
      IF (IEND.EQ.2) THEN
          CM(9) = EB + EB1 - CM(9) + ECOLL1
      END IF
*
*       Set relative energy production and total energy change.
      DB = CM(9)/EBCH0
      DE = DE + ABS(ENERGY - E0)
*
*       Provide diagnostics for exchange (membership may have changed).
      IF (NAME1 + NAME2.NE.NAMEC(I1) + NAMEC(I2).AND.EB.LT.0.0) THEN
          ISUB = ISYS(5)
          TCH = T0S(ISUB) + TIMEC
          if(rank.eq.0)
     &    WRITE (6,15)  TCH, NAME1, NAME2, NAMEC(I1), NAMEC(I2), ECC0,
     &                  ECC, A0, SEMI, EB0, EB
   15     FORMAT (' EXCHANGE    T NAM E0 E A0 A EB0 EB ',
     &                          F9.2,4I6,2F7.3,1P,2E9.1,2E10.1)
      END IF
*
*       Print final binary if relative energy increase exceeds 0.1.
      IF (DB.GT.0.1.AND.SEMI.GT.0.0) THEN 
*       Scale binding energy of relative motion by initial value.
          E1 = (ENERGY - EB - EB1)/EB0
          EB = EB/ENERGY
          EB1 = EB1/ENERGY
          if(rank.eq.0)
     &    WRITE (6,20)  NAMEC(I1), NAMEC(I2), SEMI, ECC, EB, GB, G4,
     &                  EB1, E1, ET, DB
   20     FORMAT (' CHAIN BINARY','  NAM =',2I6,'  A =',1P,E8.1,
     &            '  E =',0P,F5.2,'  EB =',F5.2,'  GB =',1P,E8.1,
     &            '  G4 =',E8.1,'  EB1 =',0P,F5.2,'  E1 =',F5.2,
     &            '  ET =',F6.3,'  DB =',F5.1)
      END IF
*
*       Include output of strong interaction from chaos case.
      IF (IEND.EQ.2.AND.PMIN.LT.2.0*SEMI.AND.DB.GT.0.1) THEN
          CALL INCLIN(XX,VV,RC,VC,ALPHA)
          if(rank.eq.0)
     &    WRITE (6,25)  NAMEC(I3), ECC, ECC1, PMIN/SEMI, RCDOT/RCP,
     &                  SEMI, A1, RCP, GB, 180.0*ALPHA/3.14
   25     FORMAT (' RECOIL:    NAM E E1 PM/A RD A A1 RP GB IN ',
     &                         I6,F8.4,F6.2,F5.1,F6.1,1P,4E10.2,0P,F7.1)
      END IF
*
*       Check optional print diagnostics of chain regularization.
      IF (KZ30.GT.1.AND.IEND.EQ.2) THEN
          TCR = MASS**2.5/ABS(2.0D0*ENERGY)**1.5
          TC = TIMEC/TCR
          EC = ENERGY/CM(8)
          if(rank.eq.0)
     &    WRITE (6,30)  I1, I2, I3, I4, RB, R13, R24, DE, TC, NSTEP1,
     &                  NREG, NPERT, DB, EC
   30     FORMAT (/,' END CHAIN  ',4I3,'  RB =',1PE8.1,'  R13 =',E8.1,
     &              '  R24 =',E8.1,'  DE =',E8.1,'  TC =',0P,F5.1,'  #',
     &                 I5,I4,I3,'  DB =',F5.2,'  EC =',F6.3)
      END IF
*
   50 RETURN
*
      END
