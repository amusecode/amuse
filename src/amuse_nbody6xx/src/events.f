      SUBROUTINE EVENTS
*
*
*       Output of mass loss or tidal capture events.
*       --------------------------------------------
*
      INCLUDE 'common6.h'
      INTEGER  NTYPE(17)
      real*8 thookf,tbgbf
      external thookf,tbgbf
*
*
*       Check counter for stellar evolution events.
      IF (NMDOT.GT.0) THEN
          DO 5 J = 1,16
              NTYPE(J) = 0
    5     CONTINUE
*
          KM = 1
          DO 10 J = 1,N
              KW = KSTAR(J) + 1
              KW = MIN(KW,16)
              KW = MAX(KW,1)
              NTYPE(KW) = NTYPE(KW) + 1
              KM = MAX(KM,KW)
   10     CONTINUE
*
          if(rank.eq.0)then
          WRITE (6,15)
   15     FORMAT (/,6X,'NMDOT   NRG  NHE  NRS  NNH  NWD  NSN  NBH  NBS',
     &               '  ZMRG  ZMHE  ZMRS  ZMNH  ZMWD  ZMSN   ZMDOT',
     &               '  NTYPE')
          WRITE (6,20)  NMDOT, NRG, NHE, NRS, NNH, NWD, NSN, NBH, NBS,
     &                  ZMRG, ZMHE, ZMRS, ZMNH, ZMWD, ZMSN, ZMDOT,
     &                  (NTYPE(J),J=1,KM)
   20     FORMAT (' #4',I9,8I5,6F6.1,F8.1,I7,I6,9I4,I5,3I4)
          end if
      END IF
*
*       Determine turnoff mass at current cluster age (cf. routine STAR).
      IF (TIME.LE.0.0D0) THEN
          TURN = BODY1*ZMBAR
      ELSE
          TPHYS = (TIME + TOFF)*TSTAR
          TURN = BODY1*ZMBAR
          TURN2 = 2.0*TURN
   25     TM = MAX(zpars(8),thookf(turn))*tbgbf(turn)
          IF (TM.GT.TPHYS) THEN
              TURN = 1.01*TURN
          ELSE
              TURN = 0.985*TURN
          END IF
          IF (ABS(TM - TPHYS).GT.1.0.AND.TURN.LT.TURN2) GO TO 25
      END IF
*
*       Check output for tidal capture, collisions or coalescence.
      IF (NDISS + NCOLL + NCOAL.GT.0.OR.EGRAV.LT.0.0D0) THEN
*       Form the net energy gain in binary interactions.
          DEGRAV = EBIN + ESUB + EBESC + EMESC + EMERGE + EGRAV - EBIN0
          ZMX = BODY1*SMU
          if(rank.eq.0)then
          WRITE (6,30)
   30     FORMAT (/,5X,'NDISS  NTIDE  NSYNC  NCOLL  NCOAL  NDD  NCIRC',
     &                 '  NROCHE  NRO  NCE  NHYP  NHYPC    EBIN ',
     &                 '  EMERGE  ECOLL  EMDOT  ECDOT  EKICK  ESESC ',
     &                 '  EBESC  EMESC  DEGRAV   EBIND  MMAX')
          WRITE (6,35)  NDISS, NTIDE, NSYNC, NCOLL, NCOAL, NDD, NCIRC,
     &                  NROCHE, NRO, NCE, NHYP, NHYPC, EBIN, EMERGE,
     &                  ECOLL, EMDOT, ECDOT, EKICK, ESESC, EBESC,
     &                  EMESC, DEGRAV, E(3), ZMX
   35     FORMAT (' #5',I8,I6,3I7,I5,I7,I8,2I5,I6,I7,3F8.3,4F7.3,F8.3,
     &                  F7.3,2F8.3,F6.1)
          end if
      END IF
*
      RETURN
*
      END
