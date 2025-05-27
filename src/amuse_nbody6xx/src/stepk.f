      SUBROUTINE STEPK(DT,DTN)
*
*
*       Block time-steps.
*       -----------------
*
      INCLUDE 'common6.h'
      DATA  ONE32 /0.03125/
*     DATA  ONE16,ONE32 /0.0625D0,0.03125/
*
*
*       Truncate any large value to first block-step entry.
      K = 1
      IF (DT.GE.1.0) GO TO 10
*
*       Compare predicted step with discrete values decreasing by 1/32.
      DT1 = DTK(1)
    1 IF (DT1.GT.DT) THEN
*         DT1 = DT1*ONE16
*         K = K + 4
          DT1 = DT1*ONE32
          K = K + 5
          GO TO 1
      END IF
*
*       Increase by 2 until predicted step is exceeded (then correct level).
    4 IF (DT1.LT.DT) THEN
          DT1 = 2.0D0*DT1
          K = K - 1
          GO TO 4
      END IF
      K = K + 1
*
   10 DTN = DTK(K)
*
      RETURN
*
      END
