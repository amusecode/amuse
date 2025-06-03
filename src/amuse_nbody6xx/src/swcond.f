      SUBROUTINE SWCOND(CHECK)
*
*
*       Check switching conditions.
*       ---------------------------
*
      INCLUDE 'commonc.h'
      COMMON/SWCALL/ NCALL
      LOGICAL CHECK
      SAVE NSW
      DATA NSW/20000/
*
*
      CHECK=.FALSE.
      NCALL=NCALL+1
*       Impose switch after every NSW steps.
      IF(NCALL.GE.NSW)THEN
      NCALL=0
      CHECK=.TRUE.
      RETURN
      END IF
*
*       Inspect the chain structure (NOTE: Inverse values 1/r are used).
      ADISTI=0.5*(N-1)/RSUM
      LRI=N-1
      DO I=1,N-2
      DO J=I+2,N
      LRI=LRI+1
*
*       Do not inspect if 1/r is small.
      IF(RINV(LRI).GT.ADISTI)THEN
      IF(J-I.GT.2)THEN
*        Check for a dangerous long loop.
*     RINVMX=MAX(RINV(I-1),RINV(I),RINV(J-1),RINV(J))
      IF(I.GT.1)THEN
      RINVMX=MAX(RINV(I-1),RINV(I))
      ELSE
      RINVMX=RINV(1)
      END IF
      RINVMX=MAX(RINVMX,RINV(J-1))
      IF(J.LT.N)RINVMX=MAX(RINVMX,RINV(J))
      IF(RINV(LRI).GT.0.5*RINVMX)THEN
      CHECK=.TRUE.
      NCALL=0
      RETURN
      END IF
      ELSE
*        See whether this is a triangle with smallest size not regularized.
      IF(RINV(LRI).GT.MAX(RINV(I),RINV(I+1)))THEN
      CHECK=.TRUE.
      NCALL=0
      RETURN
      END IF
*
      END IF
      END IF
      END DO
      END DO
*
      RETURN
      END
