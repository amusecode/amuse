      SUBROUTINE OFFSET(DTOFF)
*
*
*       Offset of global times.
*       -----------------------
*
      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
*
*
*       Update the global offset time.
    1 TOFF = TOFF + DTOFF
*
*       Reduce all individual times and epochs by offset interval.
      DO 10 I = 1,NTOT
          T0(I) = T0(I) - DTOFF
          T0R(I) = T0R(I) - DTOFF
          TEV(I) = TEV(I) - DTOFF
          TEV0(I) = TEV0(I) - DTOFF
          EPOCH(I) = EPOCH(I) - DTOFF*TSTAR
   10 CONTINUE
*
*       Set new global times.
      TIME = TIME - DTOFF
      TADJ = TADJ - DTOFF
      TNEXT = TNEXT - DTOFF
      TPREV = TPREV - DTOFF
      TBLIST = TBLIST - DTOFF
      IF (KZ(19).GT.2) THEN
          TPLOT = TPLOT - DTOFF
          TMDOT = TMDOT - DTOFF
      END IF
      DO 20 I = 1,NSUB
          T0S(I) = T0S(I) - DTOFF
          TS(I) = TS(I) - DTOFF
   20 CONTINUE
*
*       See whether more reductions are needed.
      IF (TIME.GE.DTOFF) GO TO 1
*
*       Activate control indicator for new scheduling.
      IPHASE = -1
*
      RETURN
*
      END
