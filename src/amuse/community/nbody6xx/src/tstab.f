      SUBROUTINE TSTAB(I,ECC1,SEMI1,PMIN1,YFAC,ITERM)
*
*
*       Hierarchical stability time.
*       ----------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
*
*
*       Distinguish between absolute and approximate stability.
      ITERM = 0
      IF (PMIN1.GT.YFAC*PCRIT.OR.KZ(19).EQ.0.OR.KZ(27).EQ.0) THEN
          TMDIS(NMERGE+1) = 1.0D+10
      ELSE IF (NAME(I).LT.0.OR.NAME(JCOMP).LT.0) THEN
*       Skip double hierarchy (RESET2 uses standard stability criterion).
          ITERM = 1
      ELSE IF (PMIN1.GT.0.8*YFAC*PCRIT.OR.
     &        (ECC1.GT.0.99.AND.ECC1.LT.1.0)) THEN
*       Estimate time-scale for long-lived hierarchy.
          NK = 1 + 10.0*ECC1/(1.0 - ECC1)
*       Define termination time in terms of outer periods.
          TK = TWOPI*SEMI1*SQRT(ABS(SEMI1)/(BODY(I) + BODY(JCOMP)))
*       Set disruption time in new merger location for routine KSINT.
          TMDIS(NMERGE+1) = TIME + ABS(NK*TK)
*         WRITE (6,3)  YFAC,PMIN1,PCRIT,0.8*YFAC*PCRIT,NK*TK
*   3     FORMAT (' TSTAB:    YFAC PMIN1 PCR PTEST N*TK ',F6.2,1P,4E9.2)
*       Note that large NK implies termination by other means.
          IF (ECC1.LT.1.0.AND.NK.LE.1) ITERM = 1
*       Modify PCRIT so it becomes < SEMI1*(1 - ECC1)*(1 - 2*PERT).
          IF (PCRIT.LT.PMIN1) PCRIT = 0.999*PMIN1/YFAC
      ELSE IF (ECC1.LT.1.0) THEN
*       Specify termination in all other cases (but allow hyperbolic orbit).
          ITERM = 1
      END IF
*
      RETURN
*
      END
