      SUBROUTINE TOUCH(IPAIR,I1,I2,RCOLL)
*
*
*       Collision detector for KS pairs.
*       --------------------------------
*
      INCLUDE 'common6.h'
*
*
      WRITE (6,1)  NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2),
     &             RADIUS(I1), RADIUS(I2), RCOLL, R(IPAIR), H(IPAIR)
    1 FORMAT (' KS COLL    NAM K* R* RC R H ',2I6,2I4,1P,5E10.2)
*
*       Set zero radii for binary components to avoid repeated events.
      IF (H(IPAIR).LT.0.0) THEN
          RADIUS(I1) = 0.0
          RADIUS(I2) = 0.0
      END IF
*
      NCOLL = NCOLL + 1
*
      RETURN
*
      END
