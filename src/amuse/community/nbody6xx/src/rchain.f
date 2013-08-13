      SUBROUTINE RCHAIN(IQ)
*
*
*       Regularized integration.
*       ------------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  M,MIJ,Y(25),SAVEX(12),SAVEXD(12)
      LOGICAL  SWITCH,GTYPE,GTYPE0
      COMMON/CREG/  M(4),X(12),XD(12),P(12),Q(12),TIME4,ENERGY,EPSR2,
     &              XR(9),W(9),R(6),TA(6),MIJ(6),CM(10),RMAX4,TMAX,
     &              DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/TPR/   SWITCH,GTYPE,GTYPE0
      COMMON/IND6/  IND(6)
      COMMON/CLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3,IP(4)
      COMMON/CCOLL/  QK(12),PK(12),ICALL,ICOLL,NDISS4
      EQUIVALENCE  (Y(1),P(1))
      SAVE

*
*
      IF (TIME4.GT.0.0D0) GO TO 5
      DO 1 J = 1,12
          SAVEX(J) = X(J)
          SAVEXD(J) = XD(J)
    1 CONTINUE
      TSAVE = TIME4
      TSTEP0 = TSTEP
      ISTART = 0
      NEXT = 0
      JC = 0
      S = 0.0D0
    2 ZMI1 = 100.0
      ZMI2 = 100.0
      ZMI3 = 100.0
*
      GTYPE = .FALSE.
      ISTART = ISTART + 1
      RUNAV = 0.0001*RMAX4
      EPSI = EPS
      RSTAR = SQRT(EPSR2)
*
    5 DO 50 ISTEPS = 1,10000
          GTYPE0 = GTYPE
          DSO = DS
          TIME0 = TIME4
*
*       Advance the solution one step and order the distances.
          CALL DIFSY4(25,EPSI,DS,S,Y)
          NSTEP4 = NSTEP4 + 1
          CALL RSORT(R,SWITCH,IND)
*
*       Check time transformation type.
          TRIPLE = (R(1) + R(2))*(R(2) + R(3))
          IF (TRIPLE.LT.EPSR2) THEN
              GTYPE = .TRUE.
          ELSE
              GTYPE = .FALSE.
          END IF
*
*       Update smallest moment of inertia during forward integration.
          IF (TIME4.GT.TIME0.AND.JC.EQ.0) THEN
              IMIN = IND(1)
              ZMI = R(IMIN)**2
              ZMI1 = ZMI2
              ZMI2 = ZMI3
              ZMI3 = ZMI
              RUNAV = 0.9*RUNAV + 0.1*R(IMIN)
*       Check termination during last few steps (until R(IM) > SEMI).
              IF (IQ.LT.0) GO TO 100
          END IF
*
*       Switch on search indicator during iteration or just after pericentre.
          IF (ICOLL.LT.0) ICALL = 1
          IF (R(IMIN).LT.RSTAR.AND.NSTEP4.GT.NEXT) THEN
              IF (ZMI3.GT.ZMI2.AND.ZMI2.LT.ZMI1) THEN
                  ICALL = 1
              END IF
          END IF
*
*       Delay next search a few steps to avoid the same pericentre.
          IF (ICOLL.GT.0) THEN
              NEXT = NSTEP4 + 2
              GO TO 100
          END IF
*
*       Check the termination criterion.
          IM2 = IND(5)
          IF (R(IM2).GT.RMAX4.OR.TIME4.GT.TMAX) THEN
              IF (R(IMIN).LT.0.01*RUNAV) GO TO 50
              IQ = -1
              GO TO 100
          END IF
*
*       Save old time transformation if change of type or new chain.
          IF ((GTYPE.NEQV.GTYPE0).OR.SWITCH) THEN
              TPR0 = R(1)*R(2)*R(3)
              IF (GTYPE0) TPR0 = TPR0/(R(1)*R(2) + (R(1) + R(2))*R(3))
          END IF
*
*       Check the switching condition.
          IF (SWITCH) THEN
              CALL ENDREG
              CALL NEWREG
              NREG = NREG + 1
          END IF
*
*       Modify new step after change of type or new chain.
          IF ((GTYPE.NEQV.GTYPE0).OR.SWITCH) THEN
              TPR = R(1)*R(2)*R(3)
              IF (GTYPE) TPR = TPR/(R(1)*R(2) + (R(1) + R(2))*R(3))
*       Scale regularized step by ratio of old & new time derivative.
              DS = DS*TPR0/TPR
          END IF
*
*       Check rare case of zero regularized step.
          IF (ISTEPS.LT.3.AND.DS.EQ.0.0D0) DS = 1.0E-2*DSO
          IF (DS.EQ.0.0D0) GO TO 120
   50 CONTINUE
*
  100 RETURN 
*
*       Make one restart with reduced tolerance.
  120 IF (ISTART.GT.1) RETURN
      if(rank.eq.0)
     &WRITE (6,125)  NSTEP4, TIME4, R(IND(6)), RMAX4
  125 FORMAT (5X,' RCHAIN RESTART:  # T R6 RMAX4 ',I5,3F12.6)
      TIME4 = TSAVE
      DO 130 J = 1,12
          X(J) = SAVEX(J)
          XD(J) = SAVEXD(J)
  130 CONTINUE
      SWITCH = .FALSE.
      CALL NEWREG
      EPSI = 0.1*EPS
      DS = 0.1*TSTEP0/(R(1)*R(2)*R(3))
      if(rank.eq.0)
     &WRITE (6,135)  ICALL,JC,ICOLL,DS
  135 FORMAT (' RESTART:   ICALL JC ICOLL FACM DS ',3I4,1P,E9.1)
      ICALL = 0
      ICOLL = 0
      JC = 0
      GO TO 2
*
      END
