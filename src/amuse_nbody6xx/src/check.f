      SUBROUTINE CHECK(DE)
*
*
*       Error check and restart.
*       ------------------------
*
      INCLUDE 'common6.h'
*
*
*       See whether output intervals should be increased (at main output).
      IF (KZ(32).GT.0.AND.TIME.GE.TNEXT - 20.0*DTMIN) THEN
*       Check current energy level (factor of 2) for possible increase.
          K = KZ(32)
          ECRIT = 0.25/2.0**K
          IF (ABS(E(3)).LT.ECRIT) THEN
*       Define dynamical crossing time in case energy is near zero.
              TDYN = 2.0*RSCALE/SQRT(2.0*ZKIN/ZMASS)
              IF (2.0*DTADJ.GT.TDYN.OR.TIME.LE.0.0D0) GO TO 5
              DTADJ = 2.0*DTADJ
              DELTAT = 2.0*DELTAT
              QE = SQRT(2.0)*QE
              KZ(32) = KZ(32) + 1
              if(rank.eq.0)
     &        WRITE (6,1)  DTADJ, DELTAT, QE
    1         FORMAT (/,5X,'NEW INTERVALS:   DTADJ =',F6.2,
     &                     '  DELTAT =',F6.2,'  QE =',1P,E8.1)
          END IF
      END IF
*
*       Perform automatic error check & restart (option 2).
    5 DE = ABS(DE)
      ETACOR = 1.0
*
*       Skip error check after recent mass loss (otherwise reduce KZ(19)).
      IF (KZ(19).EQ.2) THEN
          KZ(19) = KZ(19) - 1
          DE = 0.0
          GO TO 70
      END IF
*
*       Check restart for large errors (two attempts permitted).
      IF (DE.LT.5.0*QE) GO TO 30
*
*       Terminate run if no further restart is allowed.
      IF (KZ(2).LE.1.OR.NDUMP.GE.2) THEN
          if(rank.eq.0)
     &    WRITE (6,10)
   10     FORMAT (/,9X,'CALCULATIONS HALTED * * *')
*       Increase NDUMP to prevent 3rd restart (safety check in routine MAIN).
          NDUMP = 2
          IF (KZ(1).NE.0.AND.KZ(2).GE.1) CALL MYDUMP(1,1)
          call flush(6)
          STOP
      END IF
*
*       Repeat the previous interval with reduced time-step parameters.
      TCOMP = CPU
      NTEMP = NDUMP
      CALL MYDUMP(0,2)
      CPU = TCOMP
      NDUMP = NTEMP + 1
*       Control variable NDUMP used to prevent a third restart.
      ETACOR = 0.5
      ETAI = ETACOR*ETAI
      ETAR = ETACOR*ETAR
      IF (KZ(17).GT.1) ETAU = ETACOR*ETAU
      DTMIN = SQRT(ETACOR)*DTMIN
      SMIN = SQRT(ETACOR)*SMIN
      if(rank.eq.0)
     &WRITE (6,20)  TIME+TOFF, ETAI, ETAR, ETAU
   20 FORMAT (/,9X,'RESTART * * *   TIME =',F8.2,'  ETAI =',F7.3,
     &                           '  ETAR =',F7.3,'  ETAU =',F7.3)
      CALL MYDUMP(1,2)
      GO TO 50
*
*       Reset counter and check optional modification of accuracy parameters.
   30 NDUMP = 0
      IF (KZ(17).EQ.0) GO TO 50
*
      IF (DE.GT.QE) THEN
*       Continue calculation but reduce the time-step parameters.
          ETACOR = SQRT(QE/DE)
          ETAI = ETACOR*ETAI
          ETAR = ETACOR*ETAR
          IF (KZ(17).GT.1) ETAU = ETACOR*ETAU
          DTMIN = SQRT(ETACOR)*DTMIN
          SMIN = SQRT(ETACOR)*SMIN
          IF (ETACOR.LT.0.99.and.rank.eq.0)
     &    WRITE (6,40)  ETAI, ETAR, ETAU
   40     FORMAT (8X,'ETAI =',F7.3,'  ETAR =',F7.3,'  ETAU =',F7.3)
      ELSE IF (DE.LT.0.2*QE) THEN
*       Increase the time-step parameters (up to initial value only).
          IF (TIME.GT.0.0D0) THEN
              ETACOR = MIN(1.2D0,ETA0/ETAI)
              ETAI = ETACOR*ETAI
              ETAR = ETACOR*ETAR
              IF (KZ(17).GT.1) ETAU = ETACOR*ETAU
              DTMIN = SQRT(ETACOR)*DTMIN
              SMIN = SQRT(ETACOR)*SMIN
              IF (rank.eq.0.and.ETACOR.GT.1.01) 
     &        WRITE (6,40)  ETAI, ETAR, ETAU
          END IF
      END IF
*
*       See whether the time-steps should be reduced (Note: KZ(2) > 2).
   50 IF (ETACOR.LT.1.0.AND.KZ(2).GT.2) THEN
          ETACOR = SQRT(ETACOR)
          DO 60 I = IFIRST,NTOT
          IF (DMOD(T0(I),0.5D0*STEP(I)).EQ.0.0D0) THEN
              STEP(I) = 0.5*STEP(I)
              TIMENW(I) = T0(I) + STEP(I)
          END IF
          IF (DMOD(T0R(I),0.5D0*STEPR(I)).EQ.0.0D0) THEN
              STEPR(I) = 0.5*STEPR(I)
          END IF
   60     CONTINUE
*
*       Set IPHASE = -1 to ensure new time-step list in INTGRT.
          IPHASE = -1
      END IF
*
   70 RETURN
*
      END
