***
      SUBROUTINE MTURN(TURN,TPHYS,ZPARS)
*
*
*       Current MS turn-off mass.
*       --------------------------------------
*
      IMPLICIT NONE
      INTEGER  I,II,IMAX
      PARAMETER(IMAX=30)
      REAL*8 TURN,TPHYS,ZPARS(20)
      REAL*8 TM,TURN2,DM,FMID,TACC
      PARAMETER(TACC=0.001D0)
      REAL*8 THOOKF,TBGBF
      EXTERNAL THOOKF,TBGBF
*
      TURN2 = 100.D0
      TM = MAX(ZPARS(8),THOOKF(TURN2))*TBGBF(TURN2)
      IF(TM.GT.TPHYS)THEN
         TURN = TURN2
         GOTO 40
      ENDIF
*
      II = 0
 25   TM = MAX(ZPARS(8),THOOKF(TURN))*TBGBF(TURN)
      IF(TM.GT.TPHYS)THEN
         IF(TPHYS.LE.0.D0.OR.TURN.GT.98.D0) GOTO 40 
         TURN = 2.D0*TURN
         II = II + 1
         GOTO 25
      ENDIF
      TURN2 = TURN
      DM = TURN
      DO 30 , I = 1,IMAX
         DM = 0.5D0*DM
         TURN = TURN2 - DM
         TM = MAX(ZPARS(8),THOOKF(TURN))*TBGBF(TURN)
         FMID = TM - TPHYS
         IF(FMID.LT.0.0) TURN2 = TURN
         IF(DM.LT.TACC.OR.ABS(FMID).LT.1.0D-14) GOTO 40
         IF(I.EQ.IMAX)THEN
            GOTO 40
         ENDIF
 30   CONTINUE
 40   CONTINUE
*
      RETURN
*
      END
***
