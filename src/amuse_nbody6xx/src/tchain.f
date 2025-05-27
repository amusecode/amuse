      SUBROUTINE TCHAIN(ISUB,TSMIN)
*
*
*       Time interval for next chain perturber or c.m.
*       ----------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
*
*
*       Include chain c.m. in time-step search.
      NPC = LISTC(1) + 2
      LISTC(NPC) = ICH
      TSMIN = 100.0
      TIMEC = TS(ISUB)
*
*       Determine time of next chain perturber or c.m. itself.
      DO 10 L = 2,NPC
          J = LISTC(L)
          TSJ = T0(J) + STEP(J)
          TSMIN = MIN(TSMIN,TSJ - TIMEC)
   10 CONTINUE
*
*       Impose half curent step as minimum time interval.
      TSMIN = MAX(TSMIN,0.5D0*STEP(ICH))
*
*       Include safety check for large c.m. step (TPREV = TBLOCK first time).
      IF (TBLOCK.GT.TPREV) THEN
          TSMIN = MIN(TSMIN,TBLOCK - TPREV)
      ELSE
          TSMIN = MIN(TSMIN,DTMIN)
      END IF
*
*       Specify next interval unless termination.
      IF (STEPS(ISUB).GT.0.0D0) THEN
          STEPS(ISUB) = TSMIN
      END IF
*
      RETURN
*
      END
