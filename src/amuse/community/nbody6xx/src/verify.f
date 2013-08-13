      SUBROUTINE VERIFY
*
*
*       Input validation.
*       -----------------
*
      INCLUDE 'common6.h'
*
*
*       Check for unreasonable input parameters (initial & restart).
      IF (N.GE.NMAX - 2.OR.NNBMAX.GT.LMAX - 3.OR.NNBOPT.GT.NNBMAX) THEN
          if(rank.eq.0)
     &    WRITE (6,10)  N, NNBMAX, NNBOPT
   10     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   N =',I6,
     &                                 ' NNBMAX =',I4' NNBOPT =',I4)
          STOP
      END IF
*
      IF (ETAI.GT.0.08.OR.ETAR.GT.0.16) THEN
          if(rank.eq.0)
     &    WRITE (6,20)  ETAI, ETAR
   20     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   ETAI =',F6.2,
     &                                                  '  ETAR =',F6.2)
          STOP
      END IF
*
      IF (ETAU.GT.0.5.OR.GMIN.GT.0.0001.OR.GMAX.GT.0.10) THEN
          if(rank.eq.0)
     &    WRITE (6,30)  ETAU, GMIN, GMAX
   30     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   ETAU =',F6.2,
     &                                 '  GMIN =',F11.7,'  GMAX =',F7.3)
          STOP
      END IF
*
*       Also check for zero or negative values.
      IF (N.LE.0.OR.NNBMAX.LE.0.OR.ETAI.LE.0.0.OR.ETAR.LE.0.0) THEN
          if(rank.eq.0)
     &    WRITE (6,40)  N, NNBMAX, ETAI, ETAR
   40     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   N =',I5,
     &                  '  NNBMAX =',I4,'  ETAI =',F6.2,'  ETAR =',F6.2)
          STOP
      END IF
*
      IF (ETAU.LE.0.0.OR.GMIN.LE.0.0.OR.GMAX.LE.0.0) THEN
          if(rank.eq.0)
     &    WRITE (6,30)  ETAU, GMIN, GMAX
          STOP
      END IF
*
      IF (DTADJ.LE.0.0.OR.DELTAT.LE.0.0.OR.QE.LE.0.0) THEN
          if(rank.eq.0)
     &    WRITE (6,50)  DTADJ, DELTAT, QE
   50     FORMAT (/,5X,'FATAL ERROR!   BAD INPUT   DTADJ =',F6.2,
     &                                '  DELTAT =',F6.2,'  QE =',1PE9.1)
          STOP
      END IF
*
      RETURN
*
      END
